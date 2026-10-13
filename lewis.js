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
 // Aktiv "Free Draw" model: starts as 2 blank atoms + 1 bond, auto-layout (no dragging, no palette), tap to select (red),
 // RIGHT slider = bonds (neighbors 1-6, auto-adds blank atoms), LEFT slider = nonbonding electrons (or bond type for a bond),
 // Change Element / double-tap = periodic-table picker, toolbar: ? / reset / undo / FC / e-.  Extras: Check, Hint, practice mode.
 const mols=opts.mols?MOLS.filter(m=>opts.mols.includes(m.f)):MOLS;let mol=opts.mol?MOLS.find(m=>m.f==opts.mol):mols[0];
 let S,hist=[],sel=null,found=new Set(),solved=false,ref,showFC=false,showE=false,lastTap={t:0,id:null},practice=!!opts.practice,pScore=[0,0];
 const start=()=>({atoms:[{id:0,el:null,parent:-1,lp:0},{id:1,el:null,parent:0,lp:0}],ord:{1:1}});
 const reset=(keepMol)=>{S=start();hist=[];sel={t:"a",id:0};found=new Set();solved=false;ref=solveRef(mol);if(opts.preload)preloadOctet()};
 function preloadOctet(){S={atoms:mol.sk.map((x,i)=>({id:i,el:x[0],parent:x[1],lp:0})),ord:{}};S.atoms.forEach(a=>{if(a.parent>=0)S.ord[a.id]=1});let left=totalVE(mol)-2*(S.atoms.length-1);S.atoms.forEach((a,i)=>{if(i){a.lp=3;left-=6}});S.atoms[0].lp=Math.max(0,left/2);sel=null}
 const snap=()=>{hist.push(JSON.stringify(S));if(hist.length>80)hist.shift()};
 const A=id=>S.atoms.find(a=>a.id==id),kids=id=>S.atoms.filter(a=>a.parent==id),deg=a=>kids(a.id).length+(a.parent>=0?1:0);
 const bondsOf=id=>{const a=A(id),r=kids(id).map(k=>({a:id,b:k.id,o:S.ord[k.id]}));if(a.parent>=0)r.push({a:a.parent,b:id,o:S.ord[id]});return r};
 const bsum=id=>bondsOf(id).reduce((s,b)=>s+b.o,0),fcOf=a=>a.el?VAL[a.el]-2*a.lp-bsum(a.id):0,eAround=a=>2*(bsum(a.id)+a.lp);
 const need=()=>totalVE(mol),used=()=>2*Object.values(S.ord).reduce((s,o)=>s+o,0)+2*S.atoms.reduce((s,a)=>s+a.lp,0);
 const target=a=>a.el=="H"?2:a.el=="B"?6:a.el=="Be"?4:8;
 root.innerHTML=`<div class="akt akQ"><div class="akQTop"><b class="akF"></b>${opts.lockMol?"":`<select class="lm">${mols.map(m=>`<option value="${m.f}" ${m==mol?"selected":""}>${fF(m.f)}</option>`).join("")}</select><button class="btn sm sec prac">🎲 Practice mode</button><span class="chip pscore"></span>`}</div>
 <p class="mut qtxt"></p><div class="qprev" title="Tap to open the drawing tool"><div class="qpin"></div><div class="qopen">✏️ Tap to draw</div></div>
 <div class="row"><button class="btn sm sec hint">💡 Hint</button></div>
 <div class="meter" style="margin-top:6px"><i></i></div><div class="vc mut" style="font-size:13px"></div><div class="msgs"></div><div class="after"></div>
 <div class="akFS hidden"><div class="akTop"><button class="akBack" title="Back">←</button><b class="akF2"></b><button class="akDone">Done</button></div>
 <div class="akMid"><div class="slw"><input type="range" class="vsl lsl" min="0" max="8" step="2" value="0"><small class="lsll">electrons</small><b class="lslv"></b></div>
 <div class="akCan"><svg class="lwSvg" viewBox="0 0 600 430"></svg><div class="tip hidden"></div><div class="akHelp hidden"></div></div>
 <div class="slw"><input type="range" class="vsl rsl" min="1" max="6" step="1" value="1"><small>bonds</small><b class="rslv"></b></div></div>
 <div class="akBar"><button class="tb hlp" title="Help">?</button><button class="tb rst" title="Reset">🗑</button><button class="tb und" title="Undo">↶</button><button class="tb fcT" title="Formal charges">FC</button><button class="tb eT" title="Electron counts">e⁻</button><button class="tb chel" title="Change Element">Change Element</button></div></div></div>`;
 const fs=root.querySelector(".akFS"),openFS=()=>{fs.classList.remove("hidden");document.body.classList.add("akOpen")},closeFS=()=>{fs.classList.add("hidden");document.body.classList.remove("akOpen");root.querySelector(".qpin").innerHTML=root.querySelector(".lwSvg").outerHTML.replace('class="lwSvg"','class="lwPrev"')};
 root.__lb={setS:x=>{S=x;draw()},get:()=>S,select:id=>{sel={t:"a",id};draw()}};root.querySelector(".qprev").onclick=openFS;root.querySelector(".akBack").onclick=closeFS;
 const q=s=>root.querySelector(s),svg=q("svg"),lsl=q(".lsl"),rsl=q(".rsl");
 function msg(list){q(".msgs").innerHTML=list.map(([t,s])=>`<div class="msg ${t} popin">${s}</div>`).join("")}
 function layout(){const P={},LP={};const r0=S.atoms.find(a=>a.parent<0);P[r0.id]=[0,0];
  const place=(a,inAng)=>{const ks=kids(a.id),k=ks.length+a.lp+(a.parent<0?0:1),step=360/Math.max(k,1);const st=a.parent<0?(ks.length+a.lp==2?180:a.lp?270:90):inAng+180;const ang=[];for(let i=0;i<k;i++)ang.push((st+i*step)%360);
   let slots=a.parent<0?ang:ang.slice(1);if(a.parent<0&&a.lp&&ks.length+a.lp>2)slots=[...ang.slice(a.lp),...ang.slice(0,a.lp)];
   ks.forEach((c,i)=>{const t=slots[i]*Math.PI/180,L=c.el=="H"?75:100;P[c.id]=[P[a.id][0]+L*Math.cos(t),P[a.id][1]+L*Math.sin(t)];place(c,slots[i])});LP[a.id]=slots.slice(ks.length)};place(r0,0);
  const xs=Object.values(P).map(p=>p[0]),ys=Object.values(P).map(p=>p[1]);const w=Math.max(...xs)-Math.min(...xs)+120,h=Math.max(...ys)-Math.min(...ys)+120,sc=Math.min(1.3,560/w,400/h),mx=(Math.max(...xs)+Math.min(...xs))/2,my=(Math.max(...ys)+Math.min(...ys))/2;
  for(const k in P)P[k]=[300+(P[k][0]-mx)*sc,215+(P[k][1]-my)*sc];return{P,LP,sc}}
 function draw(bad=[]){const{P,LP,sc}=layout();let h="";const red="#e5484d";
  S.atoms.forEach(a=>{if(a.parent<0)return;const[x1,y1]=P[a.parent],[x2,y2]=P[a.id],dx=x2-x1,dy=y2-y1,L=Math.hypot(dx,dy)||1,nx=-dy/L*5,ny=dx/L*5,sx=dx/L*21,sy=dy/L*21,o=S.ord[a.id],isSel=sel&&sel.t=="b"&&sel.id==a.id;
   h+=`<g class="bd" data-id="${a.id}"><line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="transparent" stroke-width="26"/>`+(o==1?[0]:o==2?[-1,1]:[-1.7,0,1.7]).map(f=>`<line x1="${x1+sx+nx*f}" y1="${y1+sy+ny*f}" x2="${x2-sx+nx*f}" y2="${y2-sy+ny*f}" stroke="${isSel?red:"#3b4166"}" stroke-width="3.2" stroke-linecap="round"/>`).join("")+"</g>"});
  S.atoms.forEach(a=>{const[x,y]=P[a.id],isSel=sel&&sel.t=="a"&&sel.id==a.id,hl=bad.includes(a.id);(LP[a.id]||[]).forEach(t=>{t*=Math.PI/180;const cx=x+29*Math.cos(t),cy=y+29*Math.sin(t),px=-Math.sin(t)*6,py=Math.cos(t)*6;h+=`<circle cx="${cx+px}" cy="${cy+py}" r="3.8" fill="${isSel?red:"#7a4cff"}"/><circle cx="${cx-px}" cy="${cy-py}" r="3.8" fill="${isSel?red:"#7a4cff"}"/>`});
   h+=`<g class="at" data-id="${a.id}"><circle cx="${x}" cy="${y}" r="22" fill="${a.el?(hl?"#ffe1e1":"#fff"):"#1f2540"}" stroke="${isSel?red:hl?red:a.el?"#c9cde3":"#1f2540"}" stroke-width="${isSel?4:2.5}"/>${a.el?`<text x="${x}" y="${y+7}" text-anchor="middle" font-size="21" font-weight="800" fill="${isSel?red:"#1f2540"}">${a.el}</text>`:""}`;
   if(showFC&&a.el){const f=fcOf(a);h+=`<circle cx="${x+22}" cy="${y-20}" r="11" fill="${f>0?"#e5484d":f<0?"#3fa9ff":"#9aa0bd"}"/><text x="${x+22}" y="${y-15.5}" text-anchor="middle" font-size="12" font-weight="900" fill="#fff">${f>0?"+"+f:f<0?"−"+(-f):"0"}</text>`}
   if(showE&&a.el){const e=eAround(a),ok=a.el=="H"?e==2:(BY[a.el].p>2?e>=8:e==target(a));h+=`<rect x="${x-34}" y="${y+18}" width="30" height="16" rx="8" fill="${ok?"#18a058":"#f2b705"}"/><text x="${x-19}" y="${y+30}" text-anchor="middle" font-size="11" font-weight="900" fill="#fff">${e}e⁻</text>`}h+="</g>"});
  if(mol.q)h+=`<text x="588" y="28" font-size="17" font-weight="800" fill="#7a4cff" text-anchor="end">overall charge ${mol.q<0?(mol.q<-1?-mol.q:"")+"−":"+"}</text>`;if(solved)h=`<rect x="5" y="5" width="590" height="420" rx="14" fill="none" stroke="#18a058" stroke-width="5"/>`+h;svg.innerHTML=h;
  q(".akF").innerHTML=fF(mol.f)+(mol.n?` <span class="mut" style="font-weight:600">(${mol.n})</span>`:"");q(".akF2").innerHTML=fF(mol.f);q(".qtxt").innerHTML=`Draw the Lewis structure of ${fF(mol.f)}. Tap the box to open the drawing tool; press <b>Done</b> when finished.`;setTimeout(()=>{const p=q(".qpin");if(p)p.innerHTML=svg.outerHTML.replace('class="lwSvg"','class="lwPrev"')});q(".fcT").classList.toggle("on",showFC);q(".eT").classList.toggle("on",showE);
  // sliders reflect the selection
  if(sel&&sel.t=="a"){const a=A(sel.id);lsl.disabled=rsl.disabled=false;q(".lsll").textContent="electrons";lsl.min=0;lsl.step=2;lsl.max=Math.min(8,2*Math.max(0,6-deg(a)));lsl.value=2*a.lp;q(".lslv").textContent=2*a.lp;rsl.value=deg(a);q(".rslv").textContent=deg(a);q(".chel").disabled=false}
  else if(sel&&sel.t=="b"){lsl.disabled=false;rsl.disabled=true;q(".lsll").textContent="bond type";lsl.min=1;lsl.max=3;lsl.step=1;lsl.value=S.ord[sel.id];q(".lslv").textContent=["","single","double","triple"][S.ord[sel.id]];q(".rslv").textContent="";q(".chel").disabled=true}
  else{lsl.disabled=rsl.disabled=true;q(".lslv").textContent=q(".rslv").textContent="";q(".chel").disabled=true}
  const u=used(),nd=need(),mt=q(".meter");mt.querySelector("i").style.width=Math.min(100,u/nd*100)+"%";mt.classList.toggle("over",u>nd);mt.classList.toggle("full",u==nd);
  q(".vc").innerHTML=`Valence electrons: <b>${u}</b> / ${nd} (${Object.entries(formulaCount(mol.f)).map(([e,n])=>`${n}×${VAL[e]}`).join(" + ")}${mol.q?(mol.q<0?" + ":" − ")+Math.abs(mol.q):""})`}
 const tap=(t,id)=>{const now=Date.now();const dbl=t=="a"&&lastTap.id==id&&now-lastTap.t<380;lastTap={t:now,id:t=="a"?id:null};sel={t,id};draw();if(dbl)picker()};
 svg.addEventListener("pointerup",e=>{const g=e.target.closest(".at"),b=e.target.closest(".bd");if(g)tap("a",+g.dataset.id);else if(b)tap("b",+b.dataset.id);else{sel=null;draw()}});
 rsl.oninput=()=>{if(!sel||sel.t!="a")return;const a=A(sel.id),want=+rsl.value;let d=deg(a);if(want==d)return;snap();
  while(d<want){const id=Math.max(...S.atoms.map(x=>x.id))+1;S.atoms.push({id,el:null,parent:a.id,lp:0});S.ord[id]=1;d++}
  while(d>want){const leaf=[...kids(a.id)].reverse().sort((x,y)=>(x.el?1:0)-(y.el?1:0)).find(k=>!kids(k.id).length);if(!leaf)break;S.atoms=S.atoms.filter(x=>x!=leaf);delete S.ord[leaf.id];d--}
  if(2*a.lp>2*Math.max(0,6-deg(a)))a.lp=Math.max(0,6-deg(a));draw()};
 lsl.oninput=()=>{if(!sel)return;snap();if(sel.t=="a")A(sel.id).lp=+lsl.value/2;else S.ord[sel.id]=+lsl.value;draw()};
 function picker(){if(!sel||sel.t!="a")return;const md=$("#modal");md.style.inset="8% 12% auto";md.innerHTML=`<h3 style="margin-top:0">Choose an element</h3><div class="pk" style="display:grid;grid-template-columns:repeat(18,1fr);gap:3px"></div><p><button class="btn sec cx">Cancel</button> <button class="btn sec blank">Make blank</button></p>`;
  EL.filter(e=>e.z<=86&&e.g).forEach(e=>{const d=document.createElement("div");d.className="el";d.style.gridRow=e.row;d.style.gridColumn=e.col;d.style.background=COL[e.cat];d.innerHTML=`<b>${e.s}</b>${e.z}`;d.onclick=()=>{snap();A(sel.id).el=e.s;close();draw()};md.querySelector(".pk").appendChild(d)});
  const close=()=>{md.classList.add("hidden");md.style.inset=""};md.querySelector(".cx").onclick=close;md.querySelector(".blank").onclick=()=>{snap();A(sel.id).el=null;close();draw()};md.classList.remove("hidden")}
 q(".chel").onclick=picker;q(".und").onclick=()=>{if(hist.length){S=JSON.parse(hist.pop());sel=null;draw()}};q(".rst").onclick=()=>{snap();S=start();sel={t:"a",id:0};solved=false;q(".after").innerHTML="";msg([]);draw()};
 q(".fcT").onclick=()=>{showFC=!showFC;draw()};q(".eT").onclick=()=>{showE=!showE;draw()};
 const helpTxt="<b>How it works (same as Aktiv's Free Draw tool):</b> tap an atom or bond to select it (it turns red). The <b>right slider</b> sets how many atoms bond to the selected atom; new blank atoms appear automatically. The <b>left slider</b> sets nonbonding electrons (shown as dot pairs), or the bond type when a bond is selected. Double-tap an atom or press <b>Change Element</b> to pick an element. FC± shows the calculated formal charges and e⁻ shows each atom's electron count. 🗑 resets to 2 atoms. Check and 💡 Hint are Spectra Bench extras.";
 q(".hlp").onclick=()=>{const h=q(".akHelp");h.innerHTML=helpTxt+'<br><button class="btn sm gotit2">Close</button>';h.classList.toggle("hidden");h.querySelector(".gotit2").onclick=()=>h.classList.add("hidden")};
 if(!localStorage.sbTip){const tp=q(".tip");tp.classList.remove("hidden");tp.innerHTML=`<b>First time?</b> 1) Tap an atom (it turns red). 2) Drag the right slider to add bonds. 3) Double-tap an atom to choose its element. 4) Drag the left slider for lone pairs. <button class="btn sm gotit">Got it</button>`;tp.querySelector(".gotit").onclick=()=>{localStorage.sbTip=1;tp.classList.add("hidden")}}
 const pickMol=m=>{mol=m;reset();q(".qprev").classList.remove("good","badq");q(".after").innerHTML="";msg([]);draw();opts.onMol&&opts.onMol(mol)};
 if(q(".lm"))q(".lm").onchange=e=>pickMol(MOLS.find(m=>m.f==e.target.value));
 if(q(".prac"))q(".prac").onclick=()=>{practice=true;pScore=[0,0];const m=rnd(MOLS.filter(x=>x!=mol));q(".lm").value=m.f;pickMol(m);msg([["h","🎲 Practice mode: draw each molecule. When you get it, a new one loads."]]);q(".pscore").textContent="0 / 0"};
 const sig=list=>list.map(a=>a.el+":"+a.nb.sort().join(",")).sort().join("|");
 const refSig=()=>sig(mol.sk.map((x,i)=>({el:x[0],nb:mol.sk.filter((y,j)=>j==x[1]||y[1]==i).map(y=>y[0])})));
 const mySig=()=>sig(S.atoms.map(a=>({el:a.el,nb:bondsOf(a.id).map(b=>A(b.a==a.id?b.b:b.a).el)})));
 function skeletonText(){const cs=mol.sk.map((x,i)=>i).filter(i=>isCentral(mol.sk,i));if(cs.length==1)return`Make <b>${mol.sk[0][0]}</b> the central atom: select it and set the right slider to ${mol.sk.length-1} bonds, then give each new atom its element.`;
  return"Chain the central atoms "+mol.sk.filter((x,i)=>isCentral(mol.sk,i)).map(x=>x[0]).join("–")+", then attach "+mol.sk.map((x,i)=>x[1]>=0&&!isCentral(mol.sk,i)?`${x[0]} to ${mol.sk[x[1]][0]}`:null).filter(Boolean).join(", ")+"."}
 function audit(){const errs=[],bad=[];const blanks=S.atoms.filter(a=>!a.el);if(blanks.length)return{stage:"blank",errs:[`${blanks.length} atom${blanks.length>1?"s are":" is"} still blank (black). Double-tap each one, or use Change Element.`],bad:blanks.map(a=>a.id)};
  const cnt=formulaCount(mol.f),have={};S.atoms.forEach(a=>have[a.el]=(have[a.el]||0)+1);
  if(JSON.stringify(Object.entries(cnt).sort())!=JSON.stringify(Object.entries(have).sort()))return{stage:"atoms",errs:[`${fF(mol.f)} needs exactly ${Object.entries(cnt).map(([e,n])=>n+" "+e).join(", ")}. You have ${Object.entries(have).map(([e,n])=>n+" "+e).join(", ")}.`],bad};
  if(mySig()!=refSig())return{stage:"skeleton",errs:["The connectivity isn't right. "+skeletonText()],bad};
  const u=used(),nd=need();if(u!=nd)errs.push(u<nd?`Electron count: ${u} of ${nd} placed. Add ${nd-u} more (as lone pairs or multiple bonds).`:`Too many electrons: ${u} placed, only ${nd} available. Remove ${u-nd}.`);
  S.atoms.forEach(a=>{const e=eAround(a),p=BY[a.el].p;let m=null;if(a.el=="H"){if(e!=2)m=`H has ${e} e⁻; it needs exactly 2 (one bond, no lone pairs).`}else if(a.el=="B"||a.el=="Be"){if(e>8||(e!=target(a)&&e!=8))m=`${a.el} has ${e} e⁻ (usually ${target(a)}).`}
   else if(p<=2){if(e!=8)m=`${a.el} has ${e} e⁻. Period-2 atoms need exactly 8: ${e<8?"add electrons or raise a bond order":"remove electrons"}.`}else if(e<8)m=`${a.el} has only ${e} e⁻ (needs at least 8).`;if(m){errs.push(m);bad.push(a.id)}});
  if(errs.length)return{stage:"octet",errs,bad};
  const tot=S.atoms.reduce((s,a)=>s+Math.abs(fcOf(a)),0);if(ref&&tot>ref.tot)return{stage:"minfc",errs:[`Valid, but Σ|FC| = ${tot}, and ${ref.tot} is possible. ${BY[mol.sk[0][0]].p>2?"The central atom is in period 3, so it can expand its octet: turn O⁻ lone pairs into double bonds.":"Move electrons into multiple bonds."} (Turn on FC± to see them.)`],bad:S.atoms.filter(a=>fcOf(a)).map(a=>a.id)};
  const sec=+S.atoms.reduce((s,a)=>s+Math.max(0,-fcOf(a))*(4-BY[a.el].en),0).toFixed(3);if(ref&&tot==ref.tot&&sec>ref.sec)return{stage:"minfc",errs:["Valid, but a better structure puts the negative formal charge on the more electronegative atom."],bad:S.atoms.filter(a=>fcOf(a)<0).map(a=>a.id)};
  return{stage:"ok",errs:[],bad}}
 q(".hint").onclick=()=>{const Au=audit();const tips={blank:"Give every atom an element: double-tap it (or Change Element). Start with the central atom, "+mol.sk[0][0]+".",atoms:skeletonText(),skeleton:skeletonText(),
  octet:used()<need()?`Select each outer atom and raise its electrons (left slider) until it has 8 (H gets none). Put leftovers on the central atom. ${need()-used()} e⁻ left.`:used()>need()?"Too many electrons: lower the left slider on the central atom first.":"The count is right, but an atom is short. Lower an outer atom's electrons by 2, then select its bond and set the type to double.",
  minfc:"Turn on FC± to see the charges. Lower them by turning lone pairs into multiple bonds.",ok:"Looks right. Press Check (or Done)!"};msg([["h","💡 "+tips[Au.stage]]]);draw(Au.bad)};
 const check=()=>{const Au=audit();if(Au.stage!="ok"){msg(Au.errs.map(e=>["e","✘ "+e]));draw(Au.bad);svg.classList.remove("shake");void svg.offsetWidth;svg.classList.add("shake");if(practice){pScore[1]++;q(".pscore").textContent=`${pScore[0]} / ${pScore[1]}`}opts.onAttempt&&opts.onAttempt(false);return false}
  const key=S.atoms.filter(a=>a.parent>=0).map(a=>a.parent+"-"+a.id+":"+S.ord[a.id]+":"+a.lp).sort().join(";")+"|"+S.atoms[0].lp;const isNew=!found.has(key);found.add(key);solved=true;sel=null;draw();confetti(40);scoreAnswer(true);
  const nres=ref?ref.count:1;let t=`<b>Correct!</b> ${mol.note||""}`;if(nres>1)t+=`<br>Resonance: ${found.size} of ${nres} equivalent structures found.${found.size<nres?" To draw another, move a double bond (and its lone pairs), then Check.":" All found! The real species is their average."}${!isNew?" (You already drew this one.)":""}`;
  msg([["g","✔ "+t]]);if(practice){pScore[0]++;pScore[1]++;q(".pscore").textContent=`${pScore[0]} / ${pScore[1]}`;setTimeout(()=>{const m=rnd(MOLS.filter(x=>x!=mol));q(".lm").value=m.f;pickMol(m)},1600);return true}
  opts.onSolved&&opts.onSolved(mol,found.size,nres);if(!opts.noFollow)followUp();return true};
 q(".akDone").onclick=()=>{closeFS();const ok=check();q(".qprev").classList.toggle("good",ok);q(".qprev").classList.toggle("badq",!ok);if(ok&&opts.onDone)opts.onDone()};
 function followUp(){const G=mol.geo,box=q(".after");const qs=(opts.follow||["eg","mg","ang","pol"]).filter(k=>k!="pol"||G.pol);const L={eg:["Electron geometry",EGS],mg:["Molecular geometry",MGS],ang:["Bond angle",ANGS],pol:["Polar or nonpolar?",["polar","nonpolar"]]};
  box.innerHTML=`<div class="card popin"><b>Follow-up${G.about?" (around "+G.about+")":""}</b>`+qs.map(k=>`<p>${L[k][0]}<br>${L[k][1].map(o=>`<button class="btn sm sec fq" data-k="${k}">${o}</button>`).join("")}<span data-r="${k}"></span></p>`).join("")+`<button class="btn sm sec sol">📘 Solution</button><div class="solbox"></div></div>`;
  box.querySelectorAll(".fq").forEach(b=>b.onclick=()=>{const k=b.dataset.k,ok=b.textContent==G[k];scoreAnswer(ok);box.querySelectorAll(`.fq[data-k=${k}]`).forEach(x=>{x.disabled=true;if(x.textContent==G[k])x.classList.replace("sec","ok");else if(x==b)x.classList.replace("sec","bad")});
   box.querySelector(`[data-r=${k}]`).innerHTML=`<div class="msg ${ok?"g":"e"}">${ok?"✔ Correct!":"✘ "+G[k]+"."} ${k=="pol"?(G.pol=="nonpolar"?"Symmetric shape, so the bond dipoles cancel.":"The bond dipoles don't cancel, so there is a net dipole."):k=="ang"&&G.ang.includes("<")?"Lone pairs squeeze the angle.":""}</div>`;opts.onAttempt&&opts.onAttempt(ok)});
  box.querySelector(".sol").onclick=()=>{box.querySelector(".solbox").innerHTML=solutionSVG(mol)+vseprMini(mol)}}
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
