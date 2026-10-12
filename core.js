// ---------- shared helpers, element data, practice generators ----------
const fF=f=>String(f).replace(/([A-Za-z)])(\d+)/g,(m,a,b)=>a+[...b].map(d=>"₀₁₂₃₄₅₆₇₈₉"[d]).join(""));
const $=s=>document.querySelector(s),$$=s=>[...document.querySelectorAll(s)];
const rnd=a=>a[Math.floor(Math.random()*a.length)],shuffle=a=>{a=[...a];for(let i=a.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[a[i],a[j]]=[a[j],a[i]]}return a};
let EL=[],BY={};
const SUPD="⁰¹²³⁴⁵⁶⁷⁸⁹";const sup=s=>String(s).replace(/([spdf])(\d+)/g,(m,a,b)=>a+[...b].map(c=>SUPD[c]).join(""));
const COL={Alkali:"#ffb3a7",Alkaline:"#ffd8a8",Transition:"#fff0a8",Poor:"#c8efb0",Metalloid:"#a8e6e2",Nonmetal:"#b5dcff",Noble:"#ddc2ff",Lanthanoid:"#ffc9e3",Actinoid:"#f7b7cf",Unknown:"#ddd"};
const CATNAME={Alkali:"alkali metal",Alkaline:"alkaline earth metal",Transition:"transition metal",Poor:"post-transition metal",Metalloid:"metalloid",Nonmetal:"nonmetal",Noble:"noble gas",Lanthanoid:"lanthanide",Actinoid:"actinide",Unknown:"unknown"};
function pos(z){const P=[2,10,18,36,54,86,118];let p=P.findIndex(x=>z<=x)+1,s=p==1?1:P[p-2]+1,o=z-s+1,g;
 if(p==1)g=z==1?1:18;else if(p<=3)g=o<=2?o:o+10;else if(p<=5)g=o;else{if(o<=2)g=o;else if(o<=17)return{p,g:null,row:p+3,col:o+1};else g=o-14}return{p,g,row:p,col:g}}
const TMION={Sc:"Sc³⁺",Ti:"Ti⁴⁺",V:"V³⁺",Cr:"Cr³⁺ (also Cr²⁺)",Mn:"Mn²⁺",Fe:"Fe²⁺ / Fe³⁺",Co:"Co²⁺ / Co³⁺",Ni:"Ni²⁺",Cu:"Cu⁺ / Cu²⁺",Zn:"Zn²⁺",Ag:"Ag⁺",Cd:"Cd²⁺",Au:"Au⁺ / Au³⁺",Hg:"Hg₂²⁺ (mercury(I)) / Hg²⁺",Sn:"Sn²⁺ / Sn⁴⁺",Pb:"Pb²⁺ / Pb⁴⁺",Tl:"Tl⁺",Bi:"Bi³⁺"};
function commonIon(e){if(TMION[e.s])return TMION[e.s];if(e.s=="H")return"H⁺ (H⁻ hydride)";const g=e.g;
 if(g==1)return e.s+"⁺";if(g==2)return e.s+"²⁺";if(g==13&&e.z>5)return e.s+"³⁺";if(g==15&&e.z<=15)return e.s+"³⁻";if(g==16&&e.z<=52)return e.s+"²⁻";if(g==17&&e.z<=53)return e.s+"⁻";return"—"}
const NG={He:"1s2",Ne:"1s2 2s2 2p6",Ar:"1s2 2s2 2p6 3s2 3p6",Kr:"1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6",Xe:"1s2 2s2 2p6 3s2 3p6 4s2 3d10 4p6 5s2 4d10 5p6"};
function cfgMap(str){str=String(str).replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹]/g,c=>SUPD.indexOf(c));const m={};str=str.replace(/\[(\w+)\]/,(x,g)=>" "+(NG[g]||"")+" ");
 (str.match(/\d[spdf]\d+/g)||[]).forEach(t=>{const k=t.slice(0,2);m[k]=(m[k]||0)+ +t.slice(2)});return m}
const sameCfg=(a,b)=>{const A=cfgMap(a),B=cfgMap(b),k=new Set([...Object.keys(A),...Object.keys(B)]);for(const x of k)if((A[x]||0)!=(B[x]||0))return false;return Object.keys(A).length>0};
const ORDER="1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p".split(" ");const CAP={s:2,p:6,d:10,f:14};
function aufbau(z){let r=[],n=z;for(const o of ORDER){if(n<=0)break;const c=Math.min(n,CAP[o[1]]);r.push(o+c);n-=c}return r.join(" ")}
function unpaired(cfg){const m=cfgMap(cfg);let u=0;for(const k in m){const orb=CAP[k[1]]/2,e=m[k];u+=e<=orb?e:2*orb-e}return u}
const VAL={H:1,He:2,Li:1,Be:2,B:3,C:4,N:5,O:6,F:7,Ne:8,Na:1,Mg:2,Al:3,Si:4,P:5,S:6,Cl:7,Ar:8,Br:7,I:7,Xe:8,Se:6,K:1,Ca:2};
function formulaCount(f){const c={};(f.replace(/[−⁻²³+]/g,"").match(/[A-Z][a-z]?\d*/g)||[]).forEach(t=>{const m=t.match(/([A-Z][a-z]?)(\d*)/);c[m[1]]=(c[m[1]]||0)+(+m[2]||1)});return c}
const totalVE=m=>Object.entries(formulaCount(m.f)).reduce((s,[e,n])=>s+VAL[e]*n,0)-m.q;
const loadEls=()=>fetch("properties.json").then(r=>r.json()).then(d=>{EL=d.slice(1).map(x=>{const z=+x.atomic,p=pos(z);const e={z,s:x.symbol,mass:x.weight,cat:x.series,en:x.electroneg?+x.electroneg:null,
 rad:x.radius&&x.radius.covalent?+x.radius.covalent:null,ie:x.ionize&&x.ionize["1"]?+x.ionize["1"]:null,ea:x.affinity!=null?+x.affinity:null,dens:x.density&&(x.density.stp||x.density.liquid),melt:x.melt,boil:x.boil,cfg:x.electronstring||"",...p};
 e.val=e.z==2?2:(e.g==null||(e.g>=3&&e.g<=12))?null:e.g<=2?e.g:e.g-10;e.ion=commonIon(e);return e});EL.forEach(e=>BY[e.s]=e)});
// ---------- state: XP, streaks, lesson progress ----------
const ST=Object.assign({xp:0,streak:0,best:0,les:{},cards:{},guide:{}},JSON.parse(localStorage.sb2||"{}"));
const save=()=>{localStorage.sb2=JSON.stringify(ST);hud()};
function hud(){const s=$("#streak"),x=$("#xp");if(s){s.textContent="🔥 "+ST.streak;x.textContent="⭐ "+ST.xp}}
function toast(t){const d=document.createElement("div");d.className="toast";d.innerHTML=t;document.body.appendChild(d);setTimeout(()=>d.remove(),1800)}
function confetti(n=70){const fx=$("#fx");const c=["#5b5bf7","#18a058","#f2b705","#ff6b8a","#3fc1ff"];for(let i=0;i<n;i++){const s=document.createElement("i");s.style.left=Math.random()*100+"vw";s.style.background=rnd(c);s.style.animationDelay=Math.random()*.5+"s";s.style.animationDuration=1.2+Math.random()+"s";fx.appendChild(s);setTimeout(()=>s.remove(),2600)}}
function scoreAnswer(ok){if(ok){ST.streak++;ST.xp+=10;if(ST.streak>ST.best)ST.best=ST.streak;if(ST.streak&&ST.streak%5==0){toast(`🔥 ${ST.streak} in a row!`);confetti(30)}}else ST.streak=0;save()}
// ---------- grading + question UI ----------
const norm=s=>String(s).toLowerCase().replace(/\s+/g," ").replace(/\s*\(\s*/g,"(").replace(/\s*\)/g,")").trim();
const nf=s=>String(s).replace(/\s+/g,"").replace(/[₀-₉]/g,c=>"₀₁₂₃₄₅₆₇₈₉".indexOf(c));
function grade(Q,x){if(Q.chk)return Q.chk(x);if(Q.num!=null){const v=parseFloat(String(x).replace(/\s/g,"").replace(/[×x]10\^?/i,"e").replace(/\*10\^/,"e"));return isFinite(v)&&Math.abs(v-Q.num)/Math.abs(Q.num)<.02}return norm(x)==norm(Q.a)}
function askUI(box,Q,done){box.innerHTML=`<div class="card popin"><div class="q">${Q.q}</div>${Q.svg||""}<div class="wbox"></div>
 ${Q.c?`<div class="choices">${Q.c.map(c=>`<button class="btn sec ch">${c}</button>`).join("")}</div>`:Q.widget?"":`<div><input class="ans" size="34" autocomplete="off" placeholder="${Q.ph||"Type your answer"}"> <button class="btn go">Check</button></div>`}<div class="fbw"></div></div>`;
 const card=box.firstElementChild;let answered=false;
 const fin=(x,okOverride)=>{if(answered)return;answered=true;const ok=okOverride!=null?okOverride:grade(Q,x);scoreAnswer(ok);
  card.querySelector(".fbw").innerHTML=`<div class="fb ${ok?"y":"n"} popin"><b class="${ok?"ok":"bad"}">${ok?rnd(["✔ Nice!","✔ Correct!","✔ Yes!","✔ Nailed it!"]):"✘ Not quite."}</b> ${ok?"":"Answer: <b>"+Q.a+"</b>"}<div>${Q.why||""}</div></div>`;
  if(!ok)card.classList.add("shake");card.querySelectorAll(".ch").forEach(b=>{b.disabled=true;if(norm(b.textContent)==norm(Q.a))b.classList.add("ok");else if(b.textContent==x)b.classList.add("bad")});done&&done(ok,x)};
 card.querySelectorAll(".ch").forEach(b=>b.onclick=()=>fin(b.textContent));const i=card.querySelector(".ans");if(i){setTimeout(()=>i.focus(),50);i.onkeydown=e=>{if(e.key=="Enter")fin(i.value)};card.querySelector(".go").onclick=()=>fin(i.value)}
 if(Q.widget)Q.widget(card.querySelector(".wbox"),ok=>fin(ok?"(correct diagram)":"(diagram)",ok))}
const mc=(q,a,choices,why,extra={})=>({q,a,c:shuffle([...new Set([a,...choices])]),why,...extra});
const bank=list=>()=>{const [q,a,w,why]=rnd(list);return mc(q,a,w,why)};
// ---------- generators ----------
const SP=()=>EL.filter(e=>e.z<=36&&!(e.z>=21&&e.z<=30));
const MAINS=()=>EL.filter(e=>e.z<=54&&e.g&&(e.g<=2||e.g>=13)&&e.g!=18);
const AXE=[["AX₂","linear","linear","180"],["AX₃","trigonal planar","trigonal planar","120"],["AX₂E","trigonal planar","bent","<120"],["AX₄","tetrahedral","tetrahedral","109.5"],["AX₃E","tetrahedral","trigonal pyramidal","<109.5"],["AX₂E₂","tetrahedral","bent","<109.5"],
 ["AX₅","trigonal bipyramidal","trigonal bipyramidal","90 & 120"],["AX₄E","trigonal bipyramidal","seesaw","<90 & <120"],["AX₃E₂","trigonal bipyramidal","T-shaped","<90"],["AX₂E₃","trigonal bipyramidal","linear","180"],["AX₆","octahedral","octahedral","90"],["AX₅E","octahedral","square pyramidal","<90"],["AX₄E₂","octahedral","square planar","90"]];
const EGS=["linear","trigonal planar","tetrahedral","trigonal bipyramidal","octahedral"],MGS=["linear","bent","trigonal planar","trigonal pyramidal","tetrahedral","seesaw","T-shaped","trigonal bipyramidal","square pyramidal","square planar","octahedral"],ANGS=["180","120","109.5","<120","<109.5","90","<90","90 & 120","<90 & <120"];
const pick=(arr,a,n=3)=>shuffle(arr.filter(x=>x!=a)).slice(0,n);
const h=6.626e-34,c=2.998e8,RH=2.179e-18;
const GEN={
light:{t:"Photon calculations",g(){const lam=200+Math.round(Math.random()*600),E=h*c/(lam*1e-9),nu=c/(lam*1e-9);
 return Math.random()<.5?{q:`What is the energy (J) of one photon with λ = ${lam} nm?`,num:E,ph:"e.g. 3.97e-19",a:E.toExponential(3)+" J",why:`E = hc/λ = (6.626×10⁻³⁴)(2.998×10⁸)/(${lam}×10⁻⁹ m)`}
 :{q:`What is the frequency (s⁻¹) of light with λ = ${lam} nm?`,num:nu,ph:"e.g. 6.00e14",a:nu.toExponential(3)+" s⁻¹",why:`ν = c/λ = 2.998×10⁸ / (${lam}×10⁻⁹ m)`}}},
wave:{t:"Waves & nodes",g(){const r=Math.random();if(r<.33){const n=1+Math.floor(Math.random()*5);return{q:`A standing wave on a string has ${n} half-wavelength${n>1?"s":""}. How many nodes does it have, counting both fixed ends?`,a:String(n+1),why:`n half-waves → n + 1 nodes, so ${n+1}.`}}
 if(r<.66){const[a,b]=shuffle([[700,"red"],[450,"blue"],[550,"green"],[400,"violet"],[300,"UV"],[1000,"IR"]]).slice(0,2);const hi=a[0]<b[0]?a:b;return mc(`Which photon carries more energy: ${a[1]} (${a[0]} nm) or ${b[1]} (${b[0]} nm)?`,`${hi[1]} (${hi[0]} nm)`,[`${a[1]} (${a[0]} nm)`,`${b[1]} (${b[0]} nm)`],"E = hc/λ: shorter λ means more energy.")}
 return bank([["If the wavelength of light doubles, its frequency…","halves",["doubles","stays the same","quadruples"],"c = λν with c fixed: λ × 2 → ν ÷ 2."],
  ["Which experiment shows light behaving as a particle?","the photoelectric effect",["diffraction through slits","interference fringes","refraction in a prism"],"Ejecting electrons needs photons above a threshold frequency: particle behavior."],
  ["Which shows light behaving as a wave?","interference / diffraction",["the photoelectric effect","line spectra of atoms","photon counting"],"Interference patterns are a wave property."],
  ["All EM radiation in a vacuum travels at…","2.998×10⁸ m/s",["3.00×10⁵ m/s","6.626×10⁻³⁴ m/s","depends on λ"],"c is the same for every wavelength."]])()}},
ewave:{t:"Electron as a wave",g:bank([["Why are electron energies in atoms quantized?","only whole-number standing waves fit",["electrons are tiny planets","the nucleus repels them","photons push them"],"A confined wave allows only certain wavelengths, so only certain energies."],
 ["Heisenberg's principle says…","position and momentum can't both be known exactly",["electrons orbit in circles","energy is not conserved","light is only a wave"],"Δx·Δp ≥ h/4π."],
 ["An orbital is best described as…","a region of high electron probability",["the exact path of an electron","a circular orbit","the nucleus"],"Orbitals come from ψ²: the probability density."],
 ["de Broglie's equation is…","λ = h/mv",["E = hν","c = λν","Δx·Δp ≥ h/4π"],"Every moving particle has a wavelength λ = h/mv."],
 ["A node in an orbital is where…","the probability of finding the electron is zero",["the electron is most likely","the nucleus is","two orbitals overlap"],"ψ = 0 at a node."]])},
debroglie:{t:"de Broglie wavelength",g(){const v=+(1+Math.random()*9).toFixed(1)*1e6,lam=h/(9.109e-31*v);return{q:`An electron (m = 9.109×10⁻³¹ kg) moves at ${(v/1e6).toFixed(1)}×10⁶ m/s. What is its de Broglie wavelength in m?`,num:lam,ph:"e.g. 1.2e-10",a:lam.toExponential(3)+" m",why:"λ = h/mv"}}},
bohr:{t:"H emission λ (Bohr)",g(){const f=rnd([1,2,3]),i=f+1+Math.floor(Math.random()*3);const dE=RH*(1/f**2-1/i**2),lam=h*c/dE*1e9;return{q:`What λ (in nm) does H emit when its electron drops from n = ${i} to n = ${f}?`,num:lam,ph:"nm",a:lam.toFixed(1)+" nm",why:`ΔE = 2.179×10⁻¹⁸ (1/${f}² − 1/${i}²) = ${dE.toExponential(3)} J; λ = hc/ΔE`}}},
hconcept:{t:"Line spectra concepts",g(){if(Math.random()<.5){const ts=shuffle([[3,2],[4,2],[2,1],[5,3],[6,2]]).slice(0,3);const L=([i,f])=>h*c/(RH*(1/f**2-1/i**2));const best=ts.reduce((a,b)=>L(a)>L(b)?a:b);const s=t=>`n=${t[0]} → n=${t[1]}`;
 return mc("Which emission transition gives the LONGEST wavelength?",s(best),ts.map(s),"The smallest energy gap gives the longest λ (λ = hc/ΔE).")}
 return bank([["H emits light when its electron…","drops to a lower n",["jumps to a higher n","leaves the nucleus","stays in n = 1"],"Emission: high n → low n releases a photon."],
 ["Why does H show lines instead of a rainbow?","only certain energy levels exist",["H is colored red","the prism is broken","electrons move continuously"],"Quantized levels give discrete photon energies."],
 ["Visible H lines (Balmer) all end on…","n = 2",["n = 1","n = 3","n = ∞"],"Drops to n = 1 are UV (Lyman)."]])()}},
qnN:{t:"n: shells",g(){const n=1+Math.floor(Math.random()*4),k=Math.random();if(k<.4)return{q:`How many orbitals are in shell n = ${n}?`,a:String(n*n),why:`n² = ${n*n}`};if(k<.8)return{q:`Max electrons in shell n = ${n}?`,a:String(2*n*n),why:`2n² = ${2*n*n}`};
 return mc("Which is NOT an allowed value of n?","0",["1","2","4"],"n = 1, 2, 3, … (positive integers).")}},
qnL:{t:"l: subshells",g(){const n=1+Math.floor(Math.random()*4),k=Math.random();if(k<.5)return{q:`List the allowed l values for n = ${n} (comma-separated).`,a:[...Array(n).keys()].join(", "),chk:x=>x.replace(/\s/g,"")==[...Array(n).keys()].join(","),why:`l = 0 … n−1`};
 const sub=rnd(["1s","2s","2p","3s","3p","3d","4s","4p","4d"]);if(k<.8)return mc(`What is l for a ${sub} electron?`,String("spd".indexOf(sub[1])),["0","1","2","3"],"s = 0, p = 1, d = 2");
 return mc("Which subshell does NOT exist?",rnd(["1p","2d","3f"]),["2p","3d","4s","4d"],"l must be ≤ n − 1.")}},
qnML:{t:"mₗ: orbitals",g(){const l=Math.floor(Math.random()*3);const vals=[];for(let m=-l;m<=l;m++)vals.push(m);if(Math.random()<.5)return{q:`How many mₗ values (orbitals) are there for l = ${l} (${"spd"[l]})?`,a:String(2*l+1),why:`2l + 1 = ${2*l+1}: ${vals.join(", ")}`};
 const bad=l+1;return mc(`Which mₗ is NOT allowed for a ${2+l}${"spd"[l]} electron?`,String(bad),vals.map(String),`mₗ runs −${l} … +${l}.`)}},
qn:{t:"Allowed quantum-number sets",g(){const n=1+Math.floor(Math.random()*4),l=Math.floor(Math.random()*(n+1)),ml=Math.floor(Math.random()*5)-2,ms=rnd(["+½","−½"]);const ok=l<=n-1&&Math.abs(ml)<=l;
 return{q:`Is the set (n = ${n}, l = ${l}, mₗ = ${ml}, mₛ = ${ms}) allowed?`,c:["allowed","not allowed"],a:ok?"allowed":"not allowed",why:ok?"l ≤ n−1 and |mₗ| ≤ l.":l>n-1?`l must be 0…${n-1}.`:`mₗ must be −${l}…+${l}.`}}},
qn2el:{t:"Quantum numbers → element",g(){const e=rnd(SP());const m=cfgMap(e.cfg),ks=Object.keys(m),last=ks.filter(k=>k[1]=="p").sort().pop()||ks.filter(k=>k[1]=="s").sort().pop();
 const n=+last[0],l=last[1]=="p"?1:0,k=m[last],o=2*l+1,ml=k<=o?-l+k-1:-l+k-o-1,ms=k<=o?"+½":"−½";
 return{q:`The last electron added has <b>n = ${n}, l = ${l}, mₗ = ${ml}, mₛ = ${ms}</b>. Which element is it?<div class="mut" style="font-size:13px">Convention: fill mₗ = −l…+l with +½ first, then −½.</div>`,a:e.s,ph:"symbol, e.g. Si",chk:x=>x.trim()==e.s,why:`${last}${sup(k)}: ${sup(e.cfg)} → ${e.s}`}}},
shapes:{t:"Orbital shapes",g(){if(Math.random()<.55&&window.orbitalSVG){const o=rnd(["s","px","py","pz","dxy","dxz","dyz","dx2y2","dz2"]);const a=o[0];return mc("Which type of orbital is this?",a=="s"?"s":a=="p"?"p":"d",["s","p","d"],a=="s"?"Sphere = s.":a=="p"?"Two lobes (dumbbell) = p.":"Four lobes (or a dumbbell with a ring) = d.",{svg:`<div class="vis" style="width:180px">${orbitalSVG(o,170)}</div>`})}
 return bank([["How many p orbitals are in a p subshell?","3",["1","5","7"],"l = 1 → mₗ = −1, 0, +1."],["How many d orbitals are in a d subshell?","5",["3","7","1"],"l = 2 → 5 values of mₗ."],["The shape of an s orbital is…","a sphere",["a dumbbell","a cloverleaf","a ring"],"l = 0: spherical."],
  ["A p orbital has a nodal plane…","through the nucleus",["outside the atom","nowhere","along both lobes"],"Probability is zero in the plane between the lobes."],["Which orbital looks like a dumbbell with a donut ring?","d_z²",["dₓᵧ","pₓ","s"],"d_z² is the odd one among the d orbitals."]])()}},
pauli:{t:"Pauli exclusion",g:bank([["Which orbital box violates Pauli?","[↑↑]",["[↑↓]","[↑ ]","[  ]"],"Two electrons in one orbital must have opposite spins."],["Max electrons in ONE orbital?","2",["1","6","8"],"Pauli: 2 per orbital, ↑↓."],
 ["Two electrons in the same atom can share…","n, l, and mₗ, but not mₛ",["all four quantum numbers","only mₛ","nothing"],"No two e⁻ share all four quantum numbers."],["Which 2p diagram violates Pauli?","[↑↓↑][↑ ][  ]",["[↑↓][↑ ][↑ ]","[↑ ][↑ ][↑ ]","[↑↓][↑↓][↑ ]"],"A box can't hold 3 electrons."],
 ["Which 2p diagram violates Hund (but not Pauli)?","[↑↓][  ][  ]",["[↑ ][↑ ][  ]","[↑ ][↑ ][↑ ]","[↑↓][↑ ][↑ ]"],"Pairing before every box has one electron breaks Hund."]])},
obox:{t:"Fill an orbital diagram",g(){const e=rnd(EL.filter(x=>x.z>=3&&x.z<=36).concat([BY.N,BY.O,BY.C,BY.Cr,BY.Fe]));return{q:`Fill the orbital diagram for <b>${e.s}</b> (Z = ${e.z}). Click a box to add an electron, and right-click to remove one.`,a:sup(e.cfg),why:`${sup(e.cfg)} · ${unpaired(e.cfg)} unpaired`,widget:(el,done)=>oboxWidget(el,e,done)}}},
unpaired:{t:"Unpaired electrons",g(){const e=rnd(EL.filter(x=>x.z<=36).concat([BY.Cr,BY.Cu,BY.Fe,BY.N,BY.O]));const u=unpaired(e.cfg);return{q:`How many unpaired electrons does a ground-state ${e.s} atom have?`,a:String(u),why:`${sup(e.cfg)} → Hund's rule gives ${u} unpaired.`}}},
order:{t:"Filling order",g(){const O=ORDER.slice(0,11),i=Math.floor(Math.random()*(O.length-1));if(Math.random()<.6)return mc(`Which subshell fills right after ${O[i]}?`,O[i+1],pick(O,O[i+1]),"Order: "+O.join(" "));
 const[a,b]=shuffle([...O]).slice(0,2);const first=O.indexOf(a)<O.indexOf(b)?a:b;return mc(`Which fills first: ${a} or ${b}?`,first,[a,b],"Order: "+O.join(" "))}},
config:{t:"Electron configuration",g(){const e=rnd(EL.filter(x=>x.z<=36&&x.z>2).concat([BY.Cr,BY.Cu,BY.Cr,BY.Cu]));return{q:`Write the noble-gas configuration of <b>${e.s}</b> (Z = ${e.z}).`,ph:"e.g. [Ne] 3s2 3p1",a:sup(e.cfg),chk:x=>sameCfg(x,e.cfg),why:sup(e.cfg)+(e.s=="Cr"||e.s=="Cu"?" ← exception: "+(e.s=="Cr"?"a half-filled":"a filled")+" d is more stable":"")}}},
trend:{t:"Periodic trends",g(){const k=rnd(["rad","ie","ea"]),L={rad:"the LARGER atomic radius",ie:"the HIGHER first ionization energy",ea:"the MORE exothermic electron affinity (releases more energy)"}[k];let a,b,t=0;
 do{a=rnd(MAINS());const same=Math.random()<.5;b=rnd(MAINS().filter(x=>x!=a&&(same?x.p==a.p:x.g==a.g)));t++}while((!b||a[k]==null||b[k]==null||a[k]==b[k]||(k=="ea"&&(a.ea<=0||b.ea<=0||Math.abs(a.ea-b.ea)<20)))&&t<300);
 const w=a[k]>b[k]?a:b;return mc(`Which has ${L}: ${a.s} or ${b.s}?`,w.s,[a.s,b.s],`${a.s}: ${a[k]}, ${b.s}: ${b[k]} ${k=="rad"?"pm":"kJ/mol"} (${a.p==b.p?"same period":"same group"}).`)}},
covalent:{t:"Covalent bonding",g:bank([["A covalent bond is…","a shared pair of electrons",["a transfer of electrons","attraction between ions","a lone pair"],"Covalent = sharing."],["Which pair most likely forms a covalent bond?","C and O",["Na and Cl","K and F","Ca and O"],"Nonmetal + nonmetal → covalent."],
 ["Bond length is the distance where…","the energy is lowest",["the atoms touch","repulsion is zero","the energy is highest"],"The minimum of the energy curve."],["What holds a covalent bond together?","both nuclei attracting the shared electrons",["magnetism","gravity","neutron sharing"],"Electrostatic attraction between the nuclei and the shared pair."],
 ["Which compound is ionic?","NaCl",["CO₂","H₂O","CH₄"],"Metal + nonmetal → electron transfer."]])},
bondtype:{t:"Bond type from ΔEN",g(){let a,b,d,t=0;do{a=rnd(MAINS().filter(x=>x.en));b=rnd(MAINS().filter(x=>x.en&&x!=a));d=Math.abs(a.en-b.en);t++}while((Math.abs(d-.4)<.15||Math.abs(d-1.8)<.15)&&t<300);
 const ans=d<.4?"nonpolar covalent":d<=1.8?"polar covalent":"ionic";return mc(`What kind of bond forms between ${a.s} (EN ${a.en}) and ${b.s} (EN ${b.en})?`,ans,["nonpolar covalent","polar covalent","ionic"],`ΔEN = ${d.toFixed(2)} (< 0.4 nonpolar, 0.4–1.8 polar, > 1.8 ionic).`)}},
name:{t:"Naming drills",g(){const[f,n,af=[],an=[]]=rnd(NAMES);if(Math.random()<.5)return{q:`Name this compound: <b>${f}</b>`,a:n,chk:x=>[n,...an].map(norm).includes(norm(x)),why:n};
 return{q:`Write the formula for <b>${n}</b>`,a:f,ph:"e.g. Ca(NO3)2",chk:x=>[f,...af].map(nf).includes(nf(x)),why:f+(af.length?` (or ${af.join(", ")})`:"")+" · charges balance; ( ) only if the polyatomic count > 1."}}},
poly:{t:"Polyatomic ions",g(){const p=rnd(POLY);return Math.random()<.5?mc(`What is the name of ${p[0]}?`,p[1],pick(POLY.map(x=>x[1]),p[1]),`${p[0]} = ${p[1]}`):mc(`What is the formula for ${p[1]}?`,p[0],pick(POLY.map(x=>x[0]),p[0]),`${p[1]} = ${p[0]}`)}},
dots:{t:"Lewis dot symbols",g(){const e=rnd(EL.filter(x=>x.z<=20||(x.z>=31&&x.z<=36)));return{q:`How many dots are in the Lewis symbol of ${e.s}?`,a:String(e.val),why:`${e.s} is in group ${e.g} → ${e.val} valence e⁻.`}}},
vecount:{t:"Count valence electrons",g(){const m=rnd(MOLS);const t=totalVE(m);return{q:`How many total valence electrons does ${fF(m.f)} have?`,a:String(t),why:Object.entries(formulaCount(m.f)).map(([e,n])=>`${n}×${VAL[e]} (${e})`).join(" + ")+(m.q?` + ${-m.q} for the charge`:"")+` = ${t}`}}},
bondorder:{t:"Single/double/triple",g(){const L=[["C–O bond in CO₂","2"],["C–N bond in HCN","3"],["N–N bond in N₂","3"],["C–C bond in C₂H₄","2"],["C–C bond in C₂H₂","3"],["O–O bond in H₂O₂","1"],["Si–S bond in SiS₂","2"],["C=O bond in CH₂O","2"],["O–H bond in H₂O","1"]];const[q,a]=rnd(L);return mc(`What is the bond order of the ${q}?`,a,["1","2","3"],`Draw it: ${q.split(" in ")[1]} needs a bond order of ${a} for every atom to complete its octet.`)}},
octet:{t:"Octet violators",g(){const m=rnd(MOLS.filter(x=>x.q==0));const cen=m.sk[0][0];const nb=m.sk.filter(x=>x[1]==0).length;const lp=(totalVE(m)-2*(m.sk.length-1)-m.sk.slice(1).reduce((s,x)=>s+(x[0]=="H"?0:6),0))/2;
 const cat={BF3:"fewer than 8",BeCl2:"fewer than 8",PCl5:"more than 8",SF4:"more than 8",ClF3:"more than 8",XeF2:"more than 8",SF6:"more than 8",XeF4:"more than 8"}[m.f.replace(/[₀-₉]/g,"")]||"exactly 8";
 return mc(`In ${fF(m.f)}, the central atom (${cen}) has…`,cat,["fewer than 8","exactly 8","more than 8"],m.note||"It follows the octet rule.")}},
fc:{t:"Formal charge",g(){const combos=[["O",2,2],["O",1,3],["O",3,1],["N",3,1],["N",4,0],["N",2,2],["C",4,0],["C",3,1],["F",1,3],["Cl",1,3],["S",2,2],["S",1,3],["S",6,0],["B",3,0],["B",4,0],["P",5,0],["Xe",2,3]];const[e,b,l]=rnd(combos);const f=VAL[e]-2*l-b;
 return{q:`A <b>${e}</b> atom makes <b>${b}</b> bond${b>1?"s":""} (count a double bond as 2) and has <b>${l}</b> lone pair${l==1?"":"s"}. What is its formal charge?`,a:String(f),chk:x=>parseInt(String(x).replace("−","-").replace("+",""))==f,ph:"e.g. -1, 0, +1",why:`FC = V − L − ½B = ${VAL[e]} − ${2*l} − ${b} = ${f}`}}},
resq:{t:"Resonance",g(){if(Math.random()<.5){const m=rnd(MOLS.filter(x=>x.res>1).concat(MOLS.filter(x=>["CO2","H2O","NH3"].includes(x.f))));return mc(`How many equivalent resonance structures does ${fF(m.f)} have?`,String(m.res),["1","2","3","4"],m.res>1?`The double bond can sit in ${m.res} equivalent positions.`:"Only one good structure.")}
 return bank([["The real structure of O₃ is…","an average of 2 resonance forms",["flipping between 2 forms","one O=O and one O–O","a triangle"],"Both O–O bonds are identical (bond order 1.5)."],["Resonance structures differ in…","electron placement only",["atom positions","number of electrons","overall charge"],"Atoms stay put; only electrons move."],
  ["The best resonance structure has…","formal charges closest to zero",["the most double bonds","positive charge on O","the most lone pairs"],"Minimize |FC| and put − on the more EN atom."],["The N–O bond order in NO₃⁻ is…","1⅓",["1","1.5","2"],"4 bonding pairs ÷ 3 bonds."],["The bond order in O₃ is…","1.5",["1","2","1⅓"],"3 bonding pairs ÷ 2 bonds."]])()}},
geoEG:{t:"Electron & molecular geometry",g(){const m=rnd(MOLS),G=m.geo,k=rnd(["eg","mg"]);return mc(`${fF(m.f)}${G.about?" (around "+G.about+")":""}: what is the ${k=="eg"?"electron":"molecular"} geometry?`,G[k],pick(k=="eg"?EGS:MGS,G[k]),`${fF(m.f)}: electron ${G.eg} · molecular ${G.mg}.`)}},
axe:{t:"AXₙEₘ → shape",g(){const[x,eg,mg,ang]=rnd(AXE);return mc(`A central atom with formula type ${x} (X = bonded atoms, E = lone pairs) has what molecular geometry?`,mg,pick(MGS,mg),`${x}: electron ${eg} → molecular ${mg} (${ang}°).`)}},
lpangle:{t:"Lone-pair squeeze",g(){const L=[["CH₄",109.5],["NH₃",107],["H₂O",104.5]];const[a,b]=shuffle(L).slice(0,2);const s=a[1]<b[1]?a:b;return mc(`Which has the SMALLER H–X–H bond angle: ${a[0]} or ${b[0]}?`,s[0],[a[0],b[0]],"More lone pairs squeeze harder: CH₄ 109.5° > NH₃ 107° > H₂O 104.5°.")}},
geoANG:{t:"Bond angles",g(){const m=rnd(MOLS),G=m.geo;return mc(`${fF(m.f)}${G.about?" (around "+G.about+")":""}: what is the bond angle?`,G.ang,pick(ANGS,G.ang),`${fF(m.f)}: ${G.mg}, ${G.ang}°.`)}},
geoPOL:{t:"Molecular polarity",g(){const m=rnd(MOLS.filter(x=>x.geo.pol));return mc(`Is ${fF(m.f)} polar or nonpolar?`,m.geo.pol,["polar","nonpolar"],m.geo.pol=="nonpolar"?`${m.geo.mg} and symmetric → the bond dipoles cancel.`:`${m.geo.mg}: the dipoles don't cancel.`)}}};
// ---------- popups ----------
function showPop(html,ev){const p=$("#pop");p.innerHTML='<span class="x">✕</span>'+html;p.classList.remove("hidden");const x=Math.min((ev.clientX||200)+12,innerWidth-360),y=Math.min((ev.clientY||200)+12,innerHeight-p.offsetHeight-10);
 p.style.left=Math.max(8,x)+"px";p.style.top=Math.max(8,y)+"px";p.querySelector(".x").onclick=()=>p.classList.add("hidden");p.querySelectorAll(".term").forEach(t=>t.onclick=e=>{e.stopPropagation();termPop(t.dataset.t,e)})}
document.addEventListener("click",e=>{if(!e.target.closest("#pop,.el,.term,tr[data-s]"))$("#pop").classList.add("hidden")});
function trendText(e){if(!e.g||(e.g>=3&&e.g<=12))return"A transition metal. The Ch 3 trends focus on main-group elements.";const per=EL.filter(x=>x.p==e.p&&x.g&&(x.g<=2||x.g>=13)&&x.g!=18);
 const rk=(arr,k)=>{const s=arr.filter(x=>x[k]!=null).sort((a,b)=>b[k]-a[k]);const i=s.indexOf(e);return i<0?"—":(i+1)+" of "+s.length};
 return`Group ${e.g}, period ${e.p}. `+(e.g==18?"A noble gas: very high IE.":`Within its period, ranked largest first: radius #${rk(per,"rad")} · IE₁ #${rk(per,"ie")}.`)}
function elPop(e,ev){ev&&ev.stopPropagation&&ev.stopPropagation();let note="";
 if(e.s=="Cr")note="<p class='ok'><b>Exception:</b> aufbau predicts [Ar] 4s² 3d⁴; the actual configuration is [Ar] 4s¹ 3d⁵ (a half-filled d is more stable).</p>";
 if(e.s=="Cu")note="<p class='ok'><b>Exception:</b> aufbau predicts [Ar] 4s² 3d⁹; the actual configuration is [Ar] 4s¹ 3d¹⁰ (a filled d is more stable).</p>";
 if(e.s=="Hg")note="<p>Mercury(I) is <b>Hg₂²⁺</b>, so mercury(I) chloride is Hg₂Cl₂.</p>";
 showPop(`<h3>${e.s} <span class="mut">· Z = ${e.z}</span></h3><div class="chip">${CATNAME[e.cat]||e.cat}</div> <span class="mut">${e.mass} g/mol</span>
 <p><b>Configuration:</b> ${sup(e.cfg)}${e.z<=36?` · ${unpaired(e.cfg)} unpaired e⁻`:""}</p>${note}
 <p>EN <b>${e.en??"—"}</b> · radius <b>${e.rad??"—"}</b> pm · IE₁ <b>${e.ie??"—"}</b> kJ/mol · EA <b>${e.ea??"—"}</b> kJ/mol<br>Valence e⁻: <b>${e.val??"varies"}</b> · Common ion: <b>${e.ion}</b></p><p class="mut" style="font-size:13px">${trendText(e)}</p>`,ev)}
function termPop(t,ev){showPop(`<h3>${t.replace("Hg2^2+","Hg₂²⁺").replace("O2^2-","O₂²⁻").replace("C2H3O2-","C₂H₃O₂⁻")}</h3><p>${TERMS[t]||""}</p>`,ev)}
const termify=s=>s.replace(/\[([^\]]+)\]/g,(m,t)=>TERMS[t]?`<span class="term" data-t="${t}">${t.replace("Hg2^2+","Hg₂²⁺").replace("O2^2-","O₂²⁻").replace("C2H3O2-","C₂H₃O₂⁻")}</span>`:m);
function bindTerms(root){root.querySelectorAll(".term").forEach(t=>t.onclick=e=>{e.stopPropagation();termPop(t.dataset.t,e)})}
GEN.fcxyz={t:"Formal charge on X/Y/Z",g(){const L=[["X",16,1,3],["Y",15,4,0],["Z",14,3,1],["X",17,1,3],["Y",16,3,1],["Z",15,2,2],["X",13,4,0]];const[s,g,b,l]=rnd(L);const v=g-10,f=v-2*l-b;
 return{q:`In a Lewis structure, atom <b>${s}</b> is in group ${g}. It has ${b} bond${b>1?"s":""} (bond order total) and ${l} lone pair${l==1?"":"s"}. What is the formal charge on ${s}?`,a:String(f),chk:x=>parseInt(String(x).replace("−","-").replace("+",""))==f,why:`Valence e⁻ = ${v}. FC = ${v} − ${2*l} − ${b} = ${f}`}}};
GEN.dipole2={t:"Smallest/largest dipole",g(){const pairs=[["HF","H","F"],["HCl","H","Cl"],["HBr","H","Br"],["HI","H","I"],["ClF","Cl","F"],["BrCl","Br","Cl"],["N₂","N","N"],["CO","C","O"],["IBr","I","Br"]];const ch=shuffle(pairs).slice(0,4);const d=p=>Math.abs(BY[p[1]].en-BY[p[2]].en);
 const big=Math.random()<.5;const ans=ch.reduce((a,b)=>(big?d(b)>d(a):d(b)<d(a))?b:a);return mc(`Which diatomic molecule has the ${big?"LARGEST":"SMALLEST"} dipole moment?`,ans[0],ch.map(x=>x[0]),ch.map(p=>`${p[0]} ΔEN ${d(p).toFixed(2)}`).join(" · ")+". A bigger ΔEN means a bigger dipole.")}};
GEN.sulfite={t:"Ion formulas",g:bank([["What is the formula of the sulfite ion?","SO₃²⁻",["SO₄²⁻","S²⁻","SO₃⁻"],"-ite has one fewer O than -ate: sulfate SO₄²⁻ → sulfite SO₃²⁻."],["What is the formula of the nitrite ion?","NO₂⁻",["NO₃⁻","N³⁻","NO₂²⁻"],"nitrate NO₃⁻ → nitrite NO₂⁻."],["What is the formula of the phosphite ion?","PO₃³⁻",["PO₄³⁻","P³⁻","PO₃²⁻"],"phosphate PO₄³⁻ → phosphite PO₃³⁻."],["What is the formula of the chlorite ion?","ClO₂⁻",["ClO₃⁻","ClO⁻","ClO₄⁻"],"chlorate ClO₃⁻ → chlorite ClO₂⁻."]])};
