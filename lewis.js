// Lewis builder v2. Interaction modeled on Aktiv's model (FINDINGS/DEEP-FINDINGS): element palette, bond/lone-pair/charge
// tools, a running valence total (updateTotalValenceCount), per-atom electron counts, undo history, clear, and help. Own code.
const isCentral=(sk,i)=>sk.some(x=>x[1]==i);
function solveRef(m){ // brute force over bond orders; outer atoms fill to an octet, leftovers go to central atoms; minimize Σ|FC|
 const n=m.sk.length,bonds=m.sk.map((x,i)=>[x[1],i]).filter(b=>b[0]>=0),T=totalVE(m);let best=null,sols=[];
 const rec=(i,ord)=>{if(i==bonds.length){const bs=Array(n).fill(0);bonds.forEach((b,k)=>{bs[b[0]]+=ord[k];bs[b[1]]+=ord[k]});let left=T-2*ord.reduce((a,b)=>a+b,0);if(left<0)return;
   const lp=Array(n).fill(0);for(let a=0;a<n;a++){if(isCentral(m.sk,a))continue;const el=m.sk[a][0];const need=el=="H"?2:8;const e=need-2*bs[a];if(e<0||e%2)return;lp[a]=e/2;left-=e}
   if(left<0||left%2)return;const cs=m.sk.map((x,i)=>i).filter(i=>isCentral(m.sk,i));
   // distribute leftovers to central atoms: fill period-2 centrals to 8 first, the rest goes to the last (period 3+) central
   for(const a of cs){const el=m.sk[a][0],p2=BY[el].p<=2;const want=p2?Math.max(0,8-2*bs[a]):left;const g=Math.min(left,want-(want%2));lp[a]+=g/2;left-=g}if(left)return;
   for(let a=0;a<n;a++){const el=m.sk[a][0],e=2*(bs[a]+lp[a]),p=BY[el].p;if(el=="H"&&e!=2)return;if(el!="H"&&p<=2&&el!="B"&&el!="Be"&&e!=8)return;if((el=="B"||el=="Be")&&e>8)return;if(p>2&&e<8)return}
   const fc=m.sk.map((x,a)=>VAL[x[0]]-2*lp[a]-bs[a]),tot=fc.reduce((s,f)=>s+Math.abs(f),0),sec=+fc.reduce((s,f,a)=>s+Math.max(0,-f)*(4-BY[m.sk[a][0]].en),0).toFixed(3);if(fc.reduce((s,f)=>s+f,0)!=m.q)return;
   const sol={ord:[...ord],lp,fc,tot,sec};if(!best||tot<best.tot||(tot==best.tot&&sec<best.sec)){best=sol;sols=[sol]}else if(tot==best.tot&&sec==best.sec)sols.push(sol);return}
  for(const o of[1,2,3])rec(i+1,[...ord,o])};rec(0,[]);return best&&{...best,bonds,count:sols.length}}
function treeLayout(sk,cx,cy,lpArr){const P={};const place=(i,inAng,L)=>{const kids=sk.map((x,j)=>j).filter(j=>sk[j][1]==i);const k=kids.length+(lpArr?lpArr[i]:0)+(sk[i][1]<0?0:1);const step=360/Math.max(k,1);
  const start=sk[i][1]<0?(kids.length==2&&!(lpArr&&lpArr[i])?180:kids.length==1?180:90):inAng+180;const ang=[];for(let q=0;q<k;q++)ang.push((start+q*step)%360);const slots=sk[i][1]<0?ang:ang.slice(1);
  kids.forEach((c,q)=>{const t=slots[q]*Math.PI/180,len=sk[c][0]=="H"?70:95;P[c]=[P[i][0]+len*Math.cos(t),P[i][1]+len*Math.sin(t)];place(c,slots[q])})};P[0]=[cx,cy];place(0,0);return P}
function LewisBuilder(root,opts={}){
 const mols=opts.mols?MOLS.filter(m=>opts.mols.includes(m.f)):MOLS;let mol=opts.mol?MOLS.find(m=>m.f==opts.mol):mols[0];
 let S,hist=[],mode="bond",sel=null,found=new Set(),solved=false,ref=solveRef(mol);
 const reset=()=>{S={atoms:[],bonds:[]};hist=[];sel=null;found=new Set();solved=false;ref=solveRef(mol);if(opts.preload)preloadOctet()};
 function preloadOctet(){const P=treeLayout(mol.sk,300,215);mol.sk.forEach((x,i)=>S.atoms.push({id:i,el:x[0],x:P[i][0],y:P[i][1],lp:0,fc:0}));mol.sk.forEach((x,i)=>{if(x[1]>=0)S.bonds.push({a:x[1],b:i,o:1})});
  const T=totalVE(mol);let left=T-2*S.bonds.length;S.atoms.forEach((a,i)=>{if(i>0){a.lp=3;left-=6}});S.atoms[0].lp=Math.max(0,left/2);S.atoms.forEach(a=>a.fc=fcOf(a))}
 const snap=()=>{hist.push(JSON.stringify(S));if(hist.length>80)hist.shift()};
 const bondsOf=id=>S.bonds.filter(b=>b.a==id||b.b==id),bsum=id=>bondsOf(id).reduce((s,b)=>s+b.o,0),fcOf=a=>VAL[a.el]-2*a.lp-bsum(a.id),eAround=a=>2*(bsum(a.id)+a.lp);
 const need=()=>totalVE(mol),used=()=>2*S.bonds.reduce((s,b)=>s+b.o,0)+2*S.atoms.reduce((s,a)=>s+a.lp,0);
 const target=a=>a.el=="H"?2:a.el=="B"?6:a.el=="Be"?4:8;
 root.innerHTML=`<div class="lw"><div><div class="row" style="align-items:center;justify-content:space-between">
 ${opts.lockMol?`<b style="font-size:18px">${fF(mol.f)}${mol.n?` <span class="mut">(${mol.n})</span>`:""}</b>`:`<label>Molecule <select class="lm">${mols.map(m=>`<option value="${m.f}" ${m==mol?"selected":""}>${fF(m.f)}</option>`).join("")}</select></label>`}
 <div class="tool">
  <button class="btn sm sec md" data-m="bond" title="Bond tool: click an atom, then another, to bond them. Click a bond to change its order.">🔗 Bond</button>
  <button class="btn sm sec md" data-m="lp+" title="Click an atom to add a lone pair">•• +</button><button class="btn sm sec md" data-m="lp-" title="Click an atom to remove a lone pair">•• −</button>
  <button class="btn sm sec md" data-m="fc+" title="Click an atom: formal charge +1">⊕</button><button class="btn sm sec md" data-m="fc-" title="Click an atom: formal charge −1">⊖</button>
  <button class="btn sm sec md" data-m="erase" title="Erase atoms or bonds">🗑</button>
  <button class="btn sm sec und" title="Undo">↶</button><button class="btn sm sec clr" title="Clear">✕</button><button class="btn sm sec hlp" title="Help">?</button></div></div>
 <div class="palette"></div><svg class="lwSvg" id="lwSvg" viewBox="0 0 600 430"></svg>
 <div class="mut" style="font-size:12.5px">Drag a tile onto the canvas, or click a tile to add it. With an atom selected, a click on a tile adds a new atom already bonded to it. Drag atoms to move them. Right-click an atom to remove a lone pair.</div></div>
 <div><div class="card"><b>Valence electrons</b><div class="meter"><i></i></div><div class="vc" style="margin-top:4px"></div>
  <label class="chk"><input type="checkbox" class="showE" checked> show electron count on each atom</label></div>
 <button class="btn ok chkb">✔ Check</button><button class="btn sec hint">💡 Hint</button><div class="msgs"></div><div class="after"></div></div></div>`;
 const svg=root.querySelector("svg"),q=s=>root.querySelector(s);
 const setMode=m=>{mode=m;root.querySelectorAll(".md").forEach(b=>b.classList.toggle("ok",b.dataset.m==m));root.querySelectorAll(".md").forEach(b=>b.classList.toggle("sec",b.dataset.m!=m))};setMode("bond");
 root.querySelectorAll(".md").forEach(b=>b.onclick=()=>setMode(b.dataset.m));
 q(".und").onclick=()=>{if(hist.length){S=JSON.parse(hist.pop());sel=null;draw()}};q(".clr").onclick=()=>{snap();S={atoms:[],bonds:[]};sel=null;draw();msg([])};
 q(".hlp").onclick=()=>msg([["h","<b>How to draw:</b> 1) Add every atom from the palette (central atom first). 2) 🔗 Click atom A, then atom B, to bond them; click a bond to cycle single → double → triple. 3) Use •• + / •• − for lone pairs. 4) Use ⊕/⊖ to mark formal charges. 5) Check. The bar shows how many valence electrons you've placed out of the total."]]);
 if(q(".lm"))q(".lm").onchange=e=>{mol=MOLS.find(m=>m.f==e.target.value);reset();q(".after").innerHTML="";msg([]);draw();opts.onMol&&opts.onMol(mol)};
 q(".showE").onchange=draw;
 function msg(list){q(".msgs").innerHTML=list.map(([t,s])=>`<div class="msg ${t} popin">${s}</div>`).join("")}
 const newId=()=>S.atoms.length?Math.max(...S.atoms.map(a=>a.id))+1:0;
 function addAtom(el,x,y,bondTo){snap();const id=newId();if(x==null){if(bondTo!=null){const p=S.atoms.find(a=>a.id==bondTo);const n=bondsOf(bondTo).length;const t=(n*(360/Math.max(4,n+1))+(n%2?180:0))*Math.PI/180;x=p.x+(el=="H"?70:92)*Math.cos(t);y=p.y+(el=="H"?70:92)*Math.sin(t)}
  else if(!S.atoms.length){x=300;y=215}else{const t=S.atoms.length*1.7;x=300+120*Math.cos(t);y=215+110*Math.sin(t)}}
  x=Math.max(25,Math.min(575,x));y=Math.max(25,Math.min(405,y));S.atoms.push({id,el,x,y,lp:0,fc:0});if(S.atoms.length==1)sel=id;if(bondTo!=null)S.bonds.push({a:bondTo,b:id,o:1});draw()}
 function lpAngles(a){const bd=bondsOf(a.id).map(b=>{const o=S.atoms.find(x=>x.id==(b.a==a.id?b.b:b.a));return Math.atan2(o.y-a.y,o.x-a.x)});const res=[];
  for(let k=0;k<a.lp;k++){let bestA=0,bestD=-1;for(let t=0;t<360;t+=10){const r=t*Math.PI/180;const d=Math.min(...[...bd,...res].map(x=>{let dd=Math.abs(x-r)%(2*Math.PI);return Math.min(dd,2*Math.PI-dd)}),9);if(d>bestD+1e-6){bestD=d;bestA=r}}res.push(bestA)}return res}
 function draw(highlight=[]){let h="";const showE=q(".showE").checked;
  S.bonds.forEach((b,i)=>{const A=S.atoms.find(x=>x.id==b.a),B=S.atoms.find(x=>x.id==b.b);const dx=B.x-A.x,dy=B.y-A.y,L=Math.hypot(dx,dy)||1,nx=-dy/L*5,ny=dx/L*5,sx=dx/L*19,sy=dy/L*19;
   const offs=b.o==1?[0]:b.o==2?[-1,1]:[-1.7,0,1.7];h+=`<g class="bd" data-i="${i}" style="cursor:pointer"><line x1="${A.x}" y1="${A.y}" x2="${B.x}" y2="${B.y}" stroke="transparent" stroke-width="20"/>`+offs.map(o=>`<line x1="${A.x+sx+nx*o}" y1="${A.y+sy+ny*o}" x2="${B.x-sx+nx*o}" y2="${B.y-sy+ny*o}" stroke="#3b4166" stroke-width="3" stroke-linecap="round"/>`).join("")+"</g>"});
  S.atoms.forEach(a=>{lpAngles(a).forEach(t=>{const cx=a.x+27*Math.cos(t),cy=a.y+27*Math.sin(t),px=-Math.sin(t)*6,py=Math.cos(t)*6;h+=`<circle cx="${cx+px}" cy="${cy+py}" r="3.6" fill="#7a4cff"/><circle cx="${cx-px}" cy="${cy-py}" r="3.6" fill="#7a4cff"/>`});
   const e=eAround(a),okE=a.el=="H"?e==2:(BY[a.el].p>2?e>=8:e==target(a)),hl=highlight.includes(a.id);
   h+=`<g class="at" data-id="${a.id}" style="cursor:pointer"><circle cx="${a.x}" cy="${a.y}" r="19" fill="${sel==a.id?"#5b5bf7":hl?"#ffe1e1":"#fff"}" stroke="${hl?"#e5484d":sel==a.id?"#3f3fc4":"#c9cde3"}" stroke-width="2.5"/>
   <text x="${a.x}" y="${a.y+7}" text-anchor="middle" font-size="20" font-weight="800" fill="${sel==a.id?"#fff":"#1f2540"}">${a.el}</text>`+
   (a.fc?`<circle cx="${a.x+19}" cy="${a.y-18}" r="10" fill="${a.fc>0?"#e5484d":"#3fa9ff"}"/><text x="${a.x+19}" y="${a.y-13.5}" text-anchor="middle" font-size="12" font-weight="900" fill="#fff">${a.fc>0?"+":"−"}${Math.abs(a.fc)>1?Math.abs(a.fc):""}</text>`:"")+
   (showE?`<text x="${a.x-21}" y="${a.y+30}" font-size="11" font-weight="800" fill="${okE?"#18a058":"#c76b00"}">${e}e⁻</text>`:"")+"</g>"});
  if(mol.q)h+=`<text x="585" y="30" font-size="18" font-weight="800" fill="#7a4cff" text-anchor="end">charge ${mol.q<0?(mol.q<-1?-mol.q:"")+"−":"+"}</text>`;
  if(solved)h=`<rect x="6" y="6" width="588" height="418" rx="14" fill="none" stroke="#18a058" stroke-width="5"/>`+h;
  svg.innerHTML=h||`<text x="300" y="215" text-anchor="middle" fill="#9aa0bd" font-size="16">Drag atoms here from the palette ↑</text>`;
  const cnt=formulaCount(mol.f),have={};S.atoms.forEach(a=>have[a.el]=(have[a.el]||0)+1);
  q(".palette").innerHTML=Object.entries(cnt).map(([e,n])=>`<div class="pchip ${(have[e]||0)>=n?"done":""}" data-e="${e}">${e}<small>${have[e]||0}/${n}</small></div>`).join("");
  const u=used(),nd=need(),mt=q(".meter");mt.querySelector("i").style.width=Math.min(100,u/nd*100)+"%";mt.classList.toggle("over",u>nd);mt.classList.toggle("full",u==nd);
  q(".vc").innerHTML=`<b style="font-size:20px">${u}</b> / ${nd} placed <span class="mut">(${Object.entries(cnt).map(([e,n])=>`${n}×${VAL[e]}`).join(" + ")}${mol.q?(mol.q<0?" + ":" − ")+Math.abs(mol.q):""})</span>`;
  bindPalette()}
 function bindPalette(){root.querySelectorAll(".pchip").forEach(ch=>{ch.onpointerdown=e=>{if(ch.classList.contains("done"))return;e.preventDefault();const el=ch.dataset.e;let moved=false,ghost=null;const sx=e.clientX,sy=e.clientY;
  const mv=ev=>{if(Math.hypot(ev.clientX-sx,ev.clientY-sy)>6){moved=true;if(!ghost){ghost=ch.cloneNode(true);ghost.style.cssText="position:fixed;pointer-events:none;z-index:99;opacity:.85";document.body.appendChild(ghost)}ghost.style.left=ev.clientX-29+"px";ghost.style.top=ev.clientY-29+"px"}};
  const up=ev=>{removeEventListener("pointermove",mv);removeEventListener("pointerup",up);ghost&&ghost.remove();if(moved){const r=svg.getBoundingClientRect();if(ev.clientX>=r.left&&ev.clientX<=r.right&&ev.clientY>=r.top&&ev.clientY<=r.bottom){const p=toSvg(ev);addAtom(el,p.x,p.y,null)}}else addAtom(el,null,null,sel)};
  addEventListener("pointermove",mv);addEventListener("pointerup",up)}})}
 const toSvg=ev=>{const r=svg.getBoundingClientRect();return{x:(ev.clientX-r.left)*600/r.width,y:(ev.clientY-r.top)*430/r.height}};
 svg.addEventListener("contextmenu",e=>{const g=e.target.closest(".at");if(g){e.preventDefault();const a=S.atoms.find(x=>x.id==+g.dataset.id);if(a.lp){snap();a.lp--;draw()}}});
 svg.addEventListener("pointerdown",e=>{const g=e.target.closest(".at"),bd=e.target.closest(".bd");if(e.button==2)return;
  if(g){const a=S.atoms.find(x=>x.id==+g.dataset.id);const st=toSvg(e);let moved=false,snapped=false;
   const mv=ev=>{const p=toSvg(ev);if(!moved&&Math.hypot(p.x-st.x,p.y-st.y)>4){moved=true}if(moved){if(!snapped){snap();snapped=true}a.x=Math.max(20,Math.min(580,p.x));a.y=Math.max(20,Math.min(410,p.y));draw()}};
   const up=()=>{removeEventListener("pointermove",mv);removeEventListener("pointerup",up);if(!moved)atomClick(a)};addEventListener("pointermove",mv);addEventListener("pointerup",up);return}
  if(bd){const b=S.bonds[+bd.dataset.i];snap();if(mode=="erase")S.bonds.splice(+bd.dataset.i,1);else b.o=b.o%3+1;draw();return}
  sel=null;draw()});
 function atomClick(a){if(mode=="bond"){if(sel==null||sel==a.id){sel=a.id;draw();return}
   const ex=S.bonds.find(b=>(b.a==sel&&b.b==a.id)||(b.b==sel&&b.a==a.id));snap();if(ex)ex.o=ex.o%3+1;else S.bonds.push({a:sel,b:a.id,o:1});sel=null;draw();return}
  snap();if(mode=="lp+")a.lp=Math.min(4,a.lp+1);if(mode=="lp-")a.lp=Math.max(0,a.lp-1);if(mode=="fc+")a.fc=Math.min(3,a.fc+1);if(mode=="fc-")a.fc=Math.max(-3,a.fc-1);
  if(mode=="erase"){S.atoms=S.atoms.filter(x=>x!=a);S.bonds=S.bonds.filter(b=>b.a!=a.id&&b.b!=a.id);if(sel==a.id)sel=null}draw()}
 const sig=list=>list.map(a=>a.el+":"+a.nb.sort().join(",")).sort().join("|");
 const refSig=()=>sig(mol.sk.map((x,i)=>({el:x[0],nb:mol.sk.filter((y,j)=>j==x[1]||y[1]==i).map(y=>y[0])})));
 const mySig=()=>sig(S.atoms.map(a=>({el:a.el,nb:bondsOf(a.id).map(b=>S.atoms.find(x=>x.id==(b.a==a.id?b.b:b.a)).el)})));
 function skeletonText(){const cs=mol.sk.map((x,i)=>i).filter(i=>isCentral(mol.sk,i));if(cs.length==1)return`Put <b>${mol.sk[0][0]}</b> in the center and bond every other atom to it.`;
  return"Chain: "+mol.sk.filter((x,i)=>isCentral(mol.sk,i)).map(x=>x[0]).join("–")+", with each outer atom on: "+mol.sk.map((x,i)=>x[1]>=0&&!isCentral(mol.sk,i)?`${x[0]}→${mol.sk[x[1]][0]}`:null).filter(Boolean).join(", ")+"."}
 function audit(){const errs=[],bad=[];const cnt=formulaCount(mol.f),have={};S.atoms.forEach(a=>have[a.el]=(have[a.el]||0)+1);
  if(JSON.stringify(Object.entries(cnt).sort())!=JSON.stringify(Object.entries(have).sort())){return{stage:"atoms",errs:[`You need exactly: ${Object.entries(cnt).map(([e,n])=>n+" "+e).join(", ")}. Check the palette counts.`],bad}}
  if(S.bonds.length<S.atoms.length-1)return{stage:"connect",errs:["Not every atom is connected yet. A molecule with "+S.atoms.length+" atoms needs at least "+(S.atoms.length-1)+" bonds."],bad};
  if(mySig()!=refSig())return{stage:"skeleton",errs:["The skeleton (which atom bonds to which) isn't right. Hint: "+skeletonText()],bad};
  const u=used(),nd=need();if(u!=nd)errs.push(u<nd?`Electron count: ${u} of ${nd} placed. Add ${nd-u} more (as lone pairs or multiple bonds).`:`Too many electrons: ${u} placed, but only ${nd} are available. Remove ${(u-nd)/2} pair(s).`);
  S.atoms.forEach(a=>{const e=eAround(a),p=BY[a.el].p;let m=null;if(a.el=="H"){if(e!=2)m=`H has ${e} e⁻; it needs exactly 2 (one bond, no lone pairs).`}
   else if(a.el=="B"||a.el=="Be"){if(e>8||(e!=target(a)&&e!=8))m=`${a.el} has ${e} e⁻ (it is usually happy with ${target(a)}).`}
   else if(p<=2){if(e!=8)m=`${a.el} has ${e} e⁻. Period-2 atoms need exactly 8 — ${e<8?"add a lone pair or a multiple bond":"remove electrons"}.`}else if(e<8)m=`${a.el} has only ${e} e⁻ (needs at least 8).`;if(m){errs.push(m);bad.push(a.id)}});
  if(errs.length)return{stage:"octet",errs,bad};
  const wrong=S.atoms.filter(a=>a.fc!=fcOf(a));if(wrong.length)return{stage:"fc",errs:[`Formal charge is off on: ${wrong.map(a=>a.el).join(", ")}. Use FC = V − L − ½B and mark it with ⊕/⊖.`],bad:wrong.map(a=>a.id)};
  const tot=S.atoms.reduce((s,a)=>s+Math.abs(a.fc),0);if(ref&&tot>ref.tot)return{stage:"minfc",errs:[`This is a valid structure, but Σ|FC| = ${tot}, and ${ref.tot} is possible. ${BY[mol.sk[0][0]].p>2?"The central atom is in period 3, so it can expand its octet: turn O⁻ lone pairs into double bonds.":"Try moving electrons into multiple bonds."}`],bad:S.atoms.filter(a=>a.fc).map(a=>a.id)};
  const sec=+S.atoms.reduce((s,a)=>s+Math.max(0,-a.fc)*(4-BY[a.el].en),0).toFixed(3);if(ref&&tot==ref.tot&&sec>ref.sec)return{stage:"minfc",errs:["Valid, but a better structure puts the negative formal charge on the more electronegative atom (O rather than "+mol.sk[0][0]+")."],bad:S.atoms.filter(a=>a.fc<0).map(a=>a.id)};
  return{stage:"ok",errs:[],bad}}
 q(".hint").onclick=()=>{const A=audit();const tips={atoms:"Start by adding all the atoms. The least electronegative non-H atom is usually central: "+mol.sk[0][0]+".",connect:skeletonText(),skeleton:skeletonText(),
  octet:used()<need()?`Give the outer atoms lone pairs until each has 8 (H gets none), then put any leftovers on the central atom. You have ${need()-used()} e⁻ left.`:used()>need()?"Too many electrons: remove lone pairs, starting with the central atom.":"Electrons are right, but an atom is short. Turn a lone pair on an outer atom into a bond: remove that pair (•• −), then click the bond to raise its order.",
  fc:"FC = V − L − ½B. Example: O with 3 lone pairs and 1 bond is 6 − 6 − 1 = −1.",minfc:"Lower the formal charges: a structure with more zeros (and − on the more electronegative atom) is better.",ok:"Looks good. Press Check!"};
  msg([["h","💡 "+tips[A.stage]]]);draw(A.bad)};
 q(".chkb").onclick=()=>{const A=audit();if(A.stage!="ok"){msg(A.errs.map(e=>["e","✘ "+e]));draw(A.bad);svg.classList.remove("shake");void svg.offsetWidth;svg.classList.add("shake");opts.onAttempt&&opts.onAttempt(false);return}
  const key=S.bonds.map(b=>[Math.min(b.a,b.b),Math.max(b.a,b.b),b.o].join("-")).sort().join(";");const isNew=!found.has(key);found.add(key);solved=true;draw();confetti(40);scoreAnswer(true);
  const nres=ref?ref.count:1;let txt=`<b>Correct!</b> ${mol.note||""}`;if(nres>1)txt+=`<br>Resonance: ${found.size} of ${nres} equivalent structures found.${found.size<nres?" Want another? Move a double bond (fix the lone pairs and charges), then Check.":" All found! The real species is their average."}${!isNew?" (You already drew this one.)":""}`;
  msg([["g","✔ "+txt]]);opts.onSolved&&opts.onSolved(mol,found.size,nres);if(!opts.noFollow)followUp()};
 function followUp(){const G=mol.geo,box=q(".after");const qs=(opts.follow||["eg","mg","ang","pol"]).filter(k=>k!="pol"||G.pol);const L={eg:["Electron geometry",EGS],mg:["Molecular geometry",MGS],ang:["Bond angle",ANGS],pol:["Polar or nonpolar?",["polar","nonpolar"]]};
  box.innerHTML=`<div class="card popin"><b>Next${G.about?" (around "+G.about+")":""}</b>`+qs.map(k=>`<p>${L[k][0]}<br>${L[k][1].map(o=>`<button class="btn sm sec fq" data-k="${k}">${o}</button>`).join("")}<span class="r" data-r="${k}"></span></p>`).join("")+`<button class="btn sm sec sol">📘 Solution</button><div class="solbox"></div></div>`;
  box.querySelectorAll(".fq").forEach(b=>b.onclick=()=>{const k=b.dataset.k,ok=b.textContent==G[k];scoreAnswer(ok);box.querySelectorAll(`.fq[data-k=${k}]`).forEach(x=>{x.disabled=true;if(x.textContent==G[k])x.classList.replace("sec","ok");else if(x==b)x.classList.replace("sec","bad")});
   box.querySelector(`[data-r=${k}]`).innerHTML=`<div class="msg ${ok?"g":"e"}">${ok?"✔ Correct!":"✘ "+G[k]+"."} ${k=="pol"?(G.pol=="nonpolar"?"The shape is symmetric, so the bond dipoles cancel.":"The bond dipoles don't cancel, so there is a net dipole."):k=="ang"&&G.ang.includes("<")?"Lone pairs squeeze the angle.":""}</div>`;opts.onAttempt&&opts.onAttempt(ok)});
  box.querySelector(".sol").onclick=()=>{box.querySelector(".solbox").innerHTML=solutionSVG(mol)+(window.vseprMini?vseprMini(mol):"")}}
 reset();draw();return{get mol(){return mol}}}
function solutionSVG(m){const r=solveRef(m);if(!r)return"";const P=treeLayout(m.sk,170,120,r.lp);let h="";
 r.bonds.forEach((b,k)=>{const[x1,y1]=P[b[0]],[x2,y2]=P[b[1]],dx=x2-x1,dy=y2-y1,L=Math.hypot(dx,dy),nx=-dy/L*4,ny=dx/L*4,sx=dx/L*14,sy=dy/L*14,o=r.ord[k];(o==1?[0]:o==2?[-1,1]:[-1.7,0,1.7]).forEach(f=>h+=`<line x1="${x1+sx+nx*f}" y1="${y1+sy+ny*f}" x2="${x2-sx+nx*f}" y2="${y2-sy+ny*f}" stroke="#3b4166" stroke-width="2.4"/>`)});
 m.sk.forEach((x,i)=>{const[cx,cy]=P[i];const used=r.bonds.filter(b=>b[0]==i||b[1]==i).map(b=>{const o=b[0]==i?b[1]:b[0];return Math.atan2(P[o][1]-cy,P[o][0]-cx)});const lps=[];
  for(let k=0;k<r.lp[i];k++){let best=0,bd=-1;for(let t=0;t<360;t+=10){const a=t*Math.PI/180;const d=Math.min(...[...used,...lps].map(u=>{let q=Math.abs(u-a)%(2*Math.PI);return Math.min(q,2*Math.PI-q)}),9);if(d>bd+1e-6){bd=d;best=a}}lps.push(best)}
  lps.forEach(t=>{const px=-Math.sin(t)*4.5,py=Math.cos(t)*4.5,qx=cx+20*Math.cos(t),qy=cy+20*Math.sin(t);h+=`<circle cx="${qx+px}" cy="${qy+py}" r="2.7" fill="#7a4cff"/><circle cx="${qx-px}" cy="${qy-py}" r="2.7" fill="#7a4cff"/>`});
  h+=`<text x="${cx}" y="${cy+6}" text-anchor="middle" font-size="17" font-weight="800">${x[0]}</text>`+(r.fc[i]?`<text x="${cx+13}" y="${cy-10}" font-size="12" font-weight="900" fill="${r.fc[i]>0?"#e5484d":"#3fa9ff"}">${r.fc[i]>0?"+":"−"}${Math.abs(r.fc[i])>1?Math.abs(r.fc[i]):""}</text>`:"")});
 return`<div class="vis popin"><b>Solution: ${fF(m.f)}</b>${r.count>1?` <span class="mut">(1 of ${r.count} resonance forms)</span>`:""}<svg viewBox="0 0 340 240" width="340" height="240">${h}</svg></div>`}
// Aktiv-review style tasks: question on the left, canvas on the right
const LTASKS={
lewis:[{mol:"H2O2",pre:true,p:"First count the valence electrons in H₂O₂, then draw its Lewis structure."},{mol:"CO2",p:"Draw the Lewis structure of CO₂."},{mol:"CH3Br",p:"Draw the Lewis structure of CH₃Br."},{mol:"SiS2",p:"Draw the Lewis structure of SiS₂ (Si is central)."},{mol:"HCN",p:"Draw HCN. It has a triple bond."},{mol:"C2H4",p:"Draw C₂H₄ (a C=C double bond)."}],
octet:[{mol:"XeF4",follow:["pol"],p:"Draw XeF₄, then decide: polar or nonpolar?"},{mol:"BF3",p:"Draw BF₃. Boron is fine with 6 electrons."},{mol:"SF4",p:"Draw SF₄ (S expands its octet)."},{mol:"PCl5",p:"Draw PCl₅."}],
fc:[{mol:"SO4²−",preload:true,p:"This sulfate structure follows the octet rule but has large formal charges. Fix it to minimize the formal charges."},{mol:"CO2",p:"Draw CO₂ and confirm every formal charge is 0."},{mol:"NO2−",p:"Draw NO₂⁻ and mark the formal charges."}],
res:[{mol:"N3−",p:"Draw the azide ion N₃⁻, then draw an equivalent resonance structure."},{mol:"O3",p:"Draw both resonance structures of O₃."},{mol:"NO3−",p:"Draw the resonance structures of NO₃⁻."},{mol:"CO3²−",p:"Draw a resonance structure of CO₃²⁻."}]};
