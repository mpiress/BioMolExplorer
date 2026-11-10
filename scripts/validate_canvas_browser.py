"""Check real canvas dragging, undo/redo and navigation in Chromium."""
import argparse
import json
import re
from pathlib import Path
import statistics
import time
from urllib.request import urlopen

from playwright.sync_api import sync_playwright


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port',type=int,default=8578)
    parser.add_argument('--output',type=Path,default=Path('/tmp/biomol-canvas-browser-20261110'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    base=f'http://127.0.0.1:{args.port}'
    def state():
        with urlopen(base+'/canvas-audit',timeout=10) as response:return json.load(response)
    with sync_playwright() as playwright:
        browser=playwright.chromium.launch(executable_path='/opt/google/chrome/chrome',headless=True,
            args=['--no-sandbox','--disable-dev-shm-usage'])
        page=browser.new_page(viewport=dict(width=1440,height=1000));errors=[];results=[]
        page.on('pageerror',lambda error:errors.append(str(error)))
        page.on('console',lambda message:errors.append(message.text) if message.type=='error' and
            not message.text.startswith('Failed to load resource:') else None)
        page.on('response',lambda response:errors.append(f'HTTP {response.status}: {response.url}') if
            response.status>=400 and not response.url.endswith('/favicon.ico') else None)
        try:
            for count in (10,50,100):
                page.goto(base+f'/#/{count}',wait_until='domcontentloaded')
                page.locator('flutter-view').wait_for(state='attached',timeout=30000)
                placeholder=page.locator('flt-semantics-placeholder')
                if placeholder.count():placeholder.evaluate('element => element.click()')
                page.get_by_text(f'Canvas ready: {count} blocks, {count-1} connections',exact=True).wait_for(timeout=30000)
                label=page.get_by_role('button',name=re.compile(r'RETRIEVAL\s+Block 001$'))
                label.wait_for(timeout=15000);box=label.bounding_box()
                if not box:raise AssertionError('Node header was not rendered')
                x=box['x']+box['width']*.5;y=box['y']+box['height']*.5
                page.mouse.move(x,y);page.mouse.down()
                for step in range(1,16):
                    page.mouse.move(x+step*5,y+step*2);page.wait_for_timeout(40)
                page.mouse.up();page.wait_for_timeout(250)
                dragged=state();position=dragged['positions']['Block 001']
                if not position or position['x']<=60:raise AssertionError('Node did not move')
                if dragged['connections']!=count-1:raise AssertionError('Connections were lost')
                if dragged['autosaves']!=1:raise AssertionError('Drag completion was not saved exactly once')
                for _ in range(40):
                    moved_box=label.bounding_box()
                    if moved_box and moved_box['x']>box['x']+20:break
                    page.wait_for_timeout(100)
                else:raise AssertionError('The model moved but the rendered node did not follow')
                page.keyboard.press('Control+z');page.wait_for_timeout(200)
                if state()['positions']['Block 001'] is not None:raise AssertionError('Undo did not restore the model')
                page.keyboard.press('Control+y');page.wait_for_timeout(200)
                if state()['positions']['Block 001']!=position:raise AssertionError('Redo did not restore the movement')
                for _ in range(40):
                    restored_box=label.bounding_box()
                    if restored_box and abs(restored_box['x']-moved_box['x'])<3:break
                    page.wait_for_timeout(100)
                else:raise AssertionError('Redo did not restore the rendered node position')
                page.get_by_role('button',name='Configure selection',exact=True).click()
                page.get_by_text('Configure Block 001',exact=True).wait_for(timeout=10000)
                page.get_by_role('button',name='Close',exact=True).click()
                page.get_by_role('button',name='Connect output of Block 001',exact=True).first.click()
                page.get_by_role('button',name=re.compile(r'^Input · .* · Block 003$')).first.click()
                page.wait_for_timeout(250)
                if state()['connections']!=count:raise AssertionError('Click connection was not rendered')
                page.keyboard.press('Control+z');page.wait_for_timeout(250)
                if state()['connections']!=count-1:raise AssertionError('Undo did not remove the connection')
                if count==10:
                    page.get_by_role('button',name='Duplicate block',exact=True).click()
                    page.wait_for_timeout(250)
                    if state()['blocks']!=count+1:raise AssertionError('Duplication failed')
                    page.get_by_role('button',name='Delete block',exact=True).click()
                    page.wait_for_timeout(250)
                    if state()['blocks']!=count:raise AssertionError('Removal failed')
                page.get_by_role('button',name='Zoom in',exact=True).click()
                page.get_by_role('button',name='Zoom out',exact=True).click()
                page.get_by_role('button',name='Reset view',exact=True).click()
                page.wait_for_timeout(250)
                page.screenshot(path=str(args.output/f'canvas-{count}.png'))
                results.append(dict(blocks=count,connections=dragged['connections'],
                    drag_events=len(dragged['drag_ms']),
                    python_drag_median_ms=statistics.median(dragged['drag_ms']),undo=True,redo=True,
                    configure=True,connect=True,zoom=True))
            page.set_viewport_size(dict(width=560,height=760));page.wait_for_timeout(250)
            page.screenshot(path=str(args.output/'canvas-100-narrow.png'))
            if state()['errors'] or errors:raise AssertionError(str(state()['errors']+errors))
        except Exception as error:
            errors.append(type(error).__name__+': '+str(error))
            try:
                page.screenshot(path=str(args.output/'failed.png'))
                (args.output/'failed.html').write_text(page.content())
            except Exception as capture_error:errors.append(str(capture_error))
            (args.output/'errors.json').write_text(json.dumps(errors,indent=2))
            raise
        finally:browser.close()
    report=dict(cases=results,errors=errors,narrow_viewport_rendered=True,
        scope='Real browser interactions; Python handler timings do not measure GPU FPS')
    (args.output/'results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
