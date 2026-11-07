"""Flet and protected molecular viewer routes share the same web origin."""
from fastapi import Request
from fastapi.responses import Response
import flet_web.fastapi as flet_fastapi
from biomolexplorer.pdb_view import PREFIX
from .branding import image_bytes


def create_web_app(session,store,viewers,secret_key=None):
    app=flet_fastapi.FastAPI()

    @app.get('/favicon.png')
    async def favicon():
        return Response(image_bytes('logo.png'),media_type='image/png',headers={'Cache-Control':'no-cache'})

    @app.get('/experimental-reports/{identifier}.html')
    async def experimental_report(identifier:str):
        # Capability URL: identifiers are random and previews expire on the
        # next cleanup cycle. Never accept paths or user-provided filenames.
        if len(identifier)!=32 or any(c not in '0123456789abcdef' for c in identifier):
            return Response(status_code=404)
        path=store.staging/'report-previews'/(identifier+'.html')
        if not path.is_file():return Response(status_code=404)
        return Response(path.read_bytes(),media_type='text/html',headers={'Cache-Control':'no-store','Content-Security-Policy':"default-src 'none'; img-src data:; style-src 'unsafe-inline'; script-src 'unsafe-inline'"})

    @app.get(PREFIX+'/{resource:path}')
    async def structure(request:Request,resource:str):
        import asyncio
        prefix=request.scope.get('root_path','')+PREFIX
        status,mime,data,headers=await asyncio.to_thread(viewers.response,PREFIX+'/'+resource,prefix)
        return Response(data,status_code=status,media_type=mime,headers=headers)

    app.mount('/',flet_fastapi.app(session,upload_dir=str(store.staging),assets_dir=None,
        max_upload_size=store.max_upload_bytes,secret_key=secret_key))
    return app


def run_web(session,store,viewers,host,port,open_browser=True,secret_key=None):
    import threading
    import webbrowser
    import uvicorn
    app=create_web_app(session,store,viewers,secret_key)
    timer=None
    if open_browser:
        address='127.0.0.1' if host in ('0.0.0.0','::') else host
        if ':' in address:address='['+address+']'
        timer=threading.Timer(1,webbrowser.open,args=(f'http://{address}:{port}',));timer.daemon=True;timer.start()
    try:uvicorn.run(app,host=host,port=port,ws='websockets-sansio',access_log=False)
    finally:
        if timer:timer.cancel()
