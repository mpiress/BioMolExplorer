"""Browser audit of every production block form at three viewport sizes."""
import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
from urllib.request import urlopen
from playwright.sync_api import sync_playwright
from biomolexplorer.catalog import TITLES


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='/tmp/biomol-settings-audit')
    args=parser.parse_args();output=Path(args.output);output.mkdir(parents=True,exist_ok=True)
    with socket.socket() as sock:
        sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
    log=(output/'server.log').open('w')
    server=subprocess.Popen([sys.executable,str(Path(__file__).with_name('audit_settings_ui.py')),
        '--port',str(port)],stdout=log,stderr=subprocess.STDOUT)
    results=[]
    try:
        for _ in range(80):
            try:
                with urlopen(f'http://127.0.0.1:{port}/',timeout=1) as response:
                    if response.status==200:break
            except OSError:time.sleep(.25)
        else:raise RuntimeError('Audit server unavailable; inspect server.log')
        with sync_playwright() as playwright:
            browser=playwright.chromium.launch(executable_path='/opt/google/chrome/chrome',headless=True,
                args=['--no-sandbox','--disable-dev-shm-usage'])
            page=browser.new_page();errors=[]
            page.on('pageerror',lambda e:errors.append(str(e)))
            page.on('console',lambda m:errors.append(m.text) if m.type=='error' else None)
            for width,height in ((1440,1000),(1024,800),(560,760)):
                page.set_viewport_size(dict(width=width,height=height))
                for operation in TITLES:
                    page.goto(f'http://127.0.0.1:{port}/#/{operation}',wait_until='domcontentloaded')
                    page.locator('flutter-view').wait_for(state='attached',timeout=15000)
                    placeholder=page.locator('flt-semantics-placeholder')
                    try:
                        if placeholder.count():placeholder.evaluate("element => element.click()")
                    except Exception:pass
                    try:page.get_by_text(TITLES[operation][0],exact=True).last.wait_for(timeout=15000)
                    except Exception:
                        page.screenshot(path=str(output/'failed.png'))
                        (output/'failed.html').write_text(page.content())
                        (output/'errors.json').write_text(json.dumps(errors,indent=2))
                        raise
                    page.wait_for_timeout(150)
                    page.screenshot(path=str(output/f'{operation}-{width}.png'))
                    if operation=='docking_vina':
                        page.mouse.move(width//2,height//2)
                        for position in range(5):
                            page.mouse.wheel(0,500);page.wait_for_timeout(120)
                            page.screenshot(path=str(output/f'vina-{width}-scroll-{position}.png'))
                    results.append(dict(operation=operation,width=width,height=height,rendered=True))
            browser.close()
        (output/'results.json').write_text(json.dumps(dict(forms=results,errors=errors),indent=2))
        print(json.dumps(dict(screens=len(results),errors=errors,output=str(output))))
    finally:
        server.terminate()
        try:server.wait(timeout=5)
        except subprocess.TimeoutExpired:server.kill();server.wait()
        log.close()


if __name__=='__main__':main()
