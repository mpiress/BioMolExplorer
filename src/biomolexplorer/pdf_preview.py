"""Render the first page of an existing PDF without changing the document."""
import shutil
import subprocess
import io
from pathlib import Path
from tempfile import TemporaryDirectory

from .visualizations import MAX_VIEW_BYTES


def pdf_preview(data):
    if not data.startswith(b'%PDF-') or len(data)>MAX_VIEW_BYTES:
        raise ValueError('PDF inválido ou maior que o limite de visualização.')
    executable=shutil.which('pdftoppm')
    if not executable:raise ValueError('O renderizador de PDF não está disponível.')
    with TemporaryDirectory(prefix='biomol-footprint-') as temp:
        source=Path(temp)/'footprint.pdf';source.write_bytes(data)
        output=Path(temp)/'preview'
        try:
            subprocess.run([executable,'-f','1','-l','1','-singlefile','-scale-to','2400',
                '-png',str(source),str(output)],check=True,capture_output=True,timeout=30)
            # PDF page margins otherwise consume a substantial part of the
            # popup. Crop only blank space, retaining a small border around
            # every axis, label and legend; the downloaded PDF is untouched.
            from PIL import Image,ImageChops
            with Image.open(output.with_suffix('.png')) as preview:
                image=preview.convert('RGB')
            bounds=ImageChops.difference(image,Image.new('RGB',image.size,'white')).getbbox()
            if bounds:
                left,top,right,bottom=bounds
                image=image.crop((max(0,left-24),max(0,top-24),min(image.width,right+24),min(image.height,bottom+24)))
            stream=io.BytesIO();image.save(stream,format='PNG')
            return stream.getvalue()
        except (subprocess.SubprocessError,OSError):
            raise ValueError('Não foi possível gerar a prévia do footprint.') from None
