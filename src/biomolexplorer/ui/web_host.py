"""Flet and protected molecular viewer routes share the same web origin."""
from fastapi import Request
from fastapi.responses import Response
import flet_web.fastapi as flet_fastapi
from biomolexplorer.pdb_view import PREFIX


def create_web_app(session,store,viewers,secret_key=None):
    app=flet_fastapi.FastAPI()

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
