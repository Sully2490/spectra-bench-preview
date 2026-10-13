// Progress sync without accounts: compressed localStorage.sb2 as a code / #sync= link
const SYNC={async enc(t){const s=new Blob([t]).stream().pipeThrough(new CompressionStream("deflate-raw"));const b=new Uint8Array(await new Response(s).arrayBuffer());return btoa(String.fromCharCode(...b)).replace(/\+/g,"-").replace(/\//g,"_").replace(/=+$/,"")},
 async dec(c){c=c.trim().replace(/^.*#sync=/,"").replace(/-/g,"+").replace(/_/g,"/");const b=Uint8Array.from(atob(c),x=>x.charCodeAt(0));return await new Response(new Blob([b]).stream().pipeThrough(new DecompressionStream("deflate-raw"))).text()},
 async apply(c){try{const t=await SYNC.dec(c);JSON.parse(t);localStorage.sb2=t;return true}catch(e){alert("That progress code didn't work. Check it was copied completely.");return false}},
 async open(){const code=await SYNC.enc(localStorage.sb2||"{}"),link=location.origin+location.pathname+"#sync="+code,m=document.getElementById("modal");
  m.innerHTML=`<div class="card" style="max-width:520px;margin:10vh auto"><h3>Sync progress</h3><p class="mut">Open this link on your other device (text or email it to yourself), or paste the code there under Import.</p>
  <b>Link</b><textarea readonly rows="3" style="width:100%" class="sl">${link}</textarea><button class="btn sm cpl">Copy link</button>
  <p><b>Code</b><textarea readonly rows="2" style="width:100%" class="sc">${code}</textarea><button class="btn sm cpc">Copy code</button></p>
  <hr><b>Import</b> (replaces this device's progress)<textarea rows="2" style="width:100%" class="si" placeholder="Paste a code or link"></textarea><button class="btn sm imp">Import</button> <button class="btn sm sec cls">Close</button></div>`;
  m.classList.remove("hidden");const q=s=>m.querySelector(s),cp=(t,b)=>{navigator.clipboard?.writeText(t);b.textContent="Copied ✔"};
  q(".cpl").onclick=e=>cp(link,e.target);q(".cpc").onclick=e=>cp(code,e.target);q(".cls").onclick=()=>m.classList.add("hidden");
  q(".imp").onclick=async()=>{if(await SYNC.apply(q(".si").value)){location.hash="Map";location.reload()}}}};
if(location.hash.startsWith("#sync=")){const c=location.hash;SYNC.apply(c).then(ok=>{if(ok)alert("Progress imported ✔");location.replace(location.pathname+"#Map");location.reload()})}
