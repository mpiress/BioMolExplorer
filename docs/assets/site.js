/* Documentation-only interactions; no application APIs or storage access. */
document.documentElement.classList.add('js');
const toggle=document.querySelector('#mobile-menu-toggle');
const nav=document.querySelector('#mobile-navigation');
function closeMenu(){nav.classList.remove('open');toggle.setAttribute('aria-expanded','false');}
toggle.addEventListener('click',()=>{const open=toggle.getAttribute('aria-expanded')!=='true';nav.classList.toggle('open',open);toggle.setAttribute('aria-expanded',String(open));});
nav.querySelectorAll('a').forEach(a=>a.addEventListener('click',closeMenu));
document.addEventListener('keydown',e=>{if(e.key==='Escape'){const open=toggle.getAttribute('aria-expanded')==='true';closeMenu();if(open)toggle.focus();}});
const dialog=document.querySelector('#example-dialog');
const image=document.querySelector('#example-full');
document.querySelectorAll('[data-example]').forEach(a=>a.addEventListener('click',e=>{if(!dialog.showModal)return;e.preventDefault();const source=a.querySelector('img');image.src=a.href;image.alt=source?.alt||'Workflow diagram';document.querySelector('#example-caption').textContent=image.alt;dialog.showModal();}));
dialog.querySelector('button').addEventListener('click',()=>dialog.close());
dialog.addEventListener('click',e=>{if(e.target===dialog){const r=dialog.getBoundingClientRect();if(e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom)dialog.close();}});
