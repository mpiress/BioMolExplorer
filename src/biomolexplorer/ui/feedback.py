"""Concise interface diagnostics; full exceptions remain in private stage logs."""
import re


def close_dialog(page,dialog):
    """Close the requested window even if a notification is above it."""
    dialog.open=False
    page.update()


def readable_error(error):
    message = str(error or '').strip()
    # Older jobs can contain the EBI error page in the exception message.
    message = re.split(r'<!doctype\b|<html\b', message, maxsplit=1, flags=re.I)[0].strip()
    message = re.sub(r'^(?:RuntimeError|ValueError|Exception):\s*', '', message)
    if len(message) > 500:
        message = message[:500].rstrip() + '… Consulte o log da etapa para os detalhes.'
    return message
