"""Isolated password-recovery state and test delivery for the local workspace."""
import hashlib
import logging
import os
import re
import secrets
import smtplib
import ssl
import threading
import time
from email.message import EmailMessage


class SMTPEmailSender:
    """Deliver recovery codes using explicit environment configuration."""

    def __init__(self,host,port,from_address,username=None,password=None,security='starttls',timeout=15):
        self.host=host;self.port=port;self.from_address=from_address
        self.username=username;self.password=password;self.security=security;self.timeout=timeout

    @classmethod
    def from_environment(cls):
        host=os.environ.get('BIOMOL_SMTP_HOST','').strip()
        username=os.environ.get('BIOMOL_SMTP_USERNAME','').strip() or None
        password=os.environ.get('BIOMOL_SMTP_PASSWORD') or None
        from_address=os.environ.get('BIOMOL_SMTP_FROM','').strip() or username
        security=os.environ.get('BIOMOL_SMTP_SECURITY','starttls').strip().lower()
        try:port=int(os.environ.get('BIOMOL_SMTP_PORT','465' if security=='ssl' else '587'))
        except ValueError:raise RuntimeError('BIOMOL_SMTP_PORT deve ser uma porta numérica válida.') from None
        if not host or not from_address:
            raise RuntimeError('Envio de recuperação não configurado. Defina BIOMOL_SMTP_HOST e BIOMOL_SMTP_FROM.')
        if security not in ('starttls','ssl','plain'):
            raise RuntimeError('BIOMOL_SMTP_SECURITY deve ser starttls, ssl ou plain.')
        if bool(username)!=bool(password):
            raise RuntimeError('Defina BIOMOL_SMTP_USERNAME e BIOMOL_SMTP_PASSWORD em conjunto.')
        return cls(host,port,from_address,username,password,security)

    def __call__(self,email,code):
        message=EmailMessage()
        message['Subject']='BioMolExplorer — código de recuperação / recovery code'
        message['From']=self.from_address;message['To']=email
        message.set_content(
            f'Seu código de recuperação do BioMolExplorer é: {code}\n\n'
            'Ele expira em 10 minutos. Se você não solicitou a alteração, ignore esta mensagem.\n\n'
            f'Your BioMolExplorer recovery code is: {code}\n\n'
            'It expires in 10 minutes. If you did not request this change, ignore this message.')
        client_class=smtplib.SMTP_SSL if self.security=='ssl' else smtplib.SMTP
        try:
            client=(client_class(self.host,self.port,timeout=self.timeout,context=ssl.create_default_context())
                if self.security=='ssl' else client_class(self.host,self.port,timeout=self.timeout))
            with client:
                if self.security=='starttls':client.starttls(context=ssl.create_default_context())
                if self.username:client.login(self.username,self.password)
                client.send_message(message)
        except (OSError,smtplib.SMTPException) as exc:
            logging.getLogger('biomolexplorer.auth').exception('Falha na entrega do código de recuperação por SMTP.',
                extra={'event':'auth.password_reset_delivery_failed','smtp_host':self.host,'smtp_port':self.port})
            raise RuntimeError('Não foi possível enviar o código por e-mail. Verifique a configuração SMTP e tente novamente.') from exc


class PasswordRecoveryService:
    """Issue short-lived OTPs and single-use reset grants without storing codes."""

    def __init__(self, store, lifetime=10 * 60, sender=None):
        self.store = store
        self.lifetime = lifetime
        self.delivery=os.environ.get('BIOMOL_PASSWORD_RESET_DELIVERY','smtp').strip().lower()
        self.sender = sender or self._configured_sender
        self.expose_code = sender is not None or self.delivery=='log'
        self._lock = threading.Lock()
        self._requests = {}
        self._grants = {}

    @staticmethod
    def _digest(value):
        return hashlib.sha256(value.encode()).hexdigest()

    @staticmethod
    def _email(value):
        email=(value or '').strip().lower()
        if not re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+',email):
            raise ValueError('Informe um e-mail válido.')
        return email

    @staticmethod
    def _mock_sender(email, code):
        logging.getLogger('biomolexplorer.auth').warning(
            'Código de recuperação BioMolExplorer para %s: %s',email,code,
            extra={'event':'auth.password_reset_code'})

    def _configured_sender(self,email,code):
        if self.delivery=='log':return self._mock_sender(email,code)
        if self.delivery!='smtp':
            raise RuntimeError('BIOMOL_PASSWORD_RESET_DELIVERY deve ser smtp ou log.')
        return SMTPEmailSender.from_environment()(email,code)

    def request(self, email):
        email=self._email(email)
        if not self.store.account_exists(email):
            raise ValueError('Não existe uma conta cadastrada com este e-mail.')
        code=f'{secrets.randbelow(100_000_000):08d}'
        try:self.sender(email,code)
        except Exception:
            with self._lock:
                self._requests.pop(email,None);self._grants.pop(email,None)
            raise
        with self._lock:
            self._requests[email]=(self._digest(code),time.time()+self.lifetime,0)
            self._grants.pop(email,None)
        return code if self.expose_code else None

    def validate(self, email, code):
        email=self._email(email);code=(code or '').strip()
        if not re.fullmatch(r'\d{8}',code):
            raise ValueError('Informe o código de 8 dígitos.')
        with self._lock:
            request=self._requests.get(email)
            if not request or request[1]<time.time():
                self._requests.pop(email,None)
                raise ValueError('O código expirou. Solicite um novo código.')
            digest,expires,attempts=request
            if attempts>=5:
                self._requests.pop(email,None)
                raise ValueError('Muitas tentativas incorretas. Solicite um novo código.')
            if not secrets.compare_digest(digest,self._digest(code)):
                self._requests[email]=(digest,expires,attempts+1)
                raise ValueError('Código incorreto.')
            grant=secrets.token_urlsafe(32)
            self._grants[email]=(self._digest(grant),time.time()+5*60)
            self._requests.pop(email,None)
            return grant

    def reset(self, email, grant, password):
        email=self._email(email)
        with self._lock:
            record=self._grants.get(email)
            if not record or record[1]<time.time() or not secrets.compare_digest(record[0],self._digest(grant or '')):
                self._grants.pop(email,None)
                raise ValueError('A autorização de recuperação expirou. Solicite um novo código.')
            self.store.reset_password(email,password)
            self._grants.pop(email,None)

