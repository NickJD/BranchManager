
const $=s=>document.querySelector(s);const selectedAssemblySequences={};let studioState=null, sampleRows=[], activeId='', lastInspection=null, cohortReview=null, inspectionRunning=false, batchRunning=false;
const SETTING_INFO=[
 ['min_quality','End-trimming Phred cutoff','Mott-style end-trimming threshold.'],
 ['min_read_length','Minimum usable read (bp)','Per-read minimum before it may contribute.'],
 ['min_final_length','Minimum final sequence (bp)','Applied to the assembled consensus or best read.'],
 ['min_mean_quality','Minimum mean quality','Mean quality for the retained read or final sequence.'],
 ['max_output_expected_errors','Maximum final expected errors','Consensus-level expected-error limit.'],
 ['warn_n_percent','N percent warning','Passing sequence requires review above this value.'],
 ['max_n_percent','N percent failure','Sequence fails above this value.'],
 ['min_overlap','Minimum read overlap (bp)','Minimum aligned overlap accepted for assembly.'],
 ['min_overlap_identity','Minimum overlap identity','Minimum identity fraction across the aligned overlap.']
];
function esc(x){return String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}
function pillClass(q){return q==='FAIL_QC'?'fail':q==='PASS_WITH_WARNINGS'?'warn':q==='PASS_HIGH_CONFIDENCE'?'pass':'neutral'}
async function api(path,body){let r=await fetch(path,body?{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)}:{});let j=await r.json();if(!r.ok)throw new Error(j.error||r.statusText);return j}
function setStatus(node,msg,error=false){node.textContent=msg;node.className=error?'error':'small'}
function renderSettings(){let p=studioState.policy,unit={min_quality:' Phred',min_read_length:' bp',min_final_length:' bp',min_mean_quality:' Phred',max_output_expected_errors:' errors',warn_n_percent:'%',max_n_percent:'%',min_overlap:' bp',min_overlap_identity:' fraction'};$('#settings').innerHTML=SETTING_INFO.map(([k,label,help])=>`<div class="setting"><label for="s_${k}">${label}</label><input id="s_${k}" type="number" step="${k.includes('identity')?'0.01':(k==='min_mean_quality'||k.includes('errors')||k.includes('percent')?'0.1':'1')}" value="${esc(p[k])}"><div class="explain">${help} Default: ${esc(p[k])}${unit[k]}.</div></div>`).join('')}
function renderSamples(){let q=$('#sampleSearch')?.value?.toLowerCase()||'';let samples=(studioState?.samples||[]).filter(s=>s.sequence_id.toLowerCase().includes(q));$('#samples').innerHTML=samples.map(s=>`<div class="sample ${s.sequence_id===activeId?'active':''}" data-id="${esc(s.sequence_id)}"><div><div class="sampleid">${esc(s.sequence_id)}</div><div class="samplemeta">${esc(s.processing_mode)} · ${s.read_count} read(s)</div></div><span class="count">${s.missing_files?'MISSING':'READY'}</span></div>`).join('')||'<div class="small" style="padding:10px">No isolates loaded.</div>';document.querySelectorAll('.sample').forEach(el=>el.onclick=()=>selectSample(el.dataset.id))}
function selectSample(id){
 if(inspectionRunning){setRunFeedback(`Inspection is still running for ${activeId}. Wait for it to finish before changing isolates.`);return}
 const s=(studioState.samples||[]).find(x=>x.sequence_id===id);if(!s)return;
 activeId=id;
 $('#sampleTaxonomyStatus').textContent=studioState.onboarding_taxonomy?.loaded?'Assistant pre-screened assignments loaded':'NOT_SCREENED — no cached assignments';
 // Processing mode is an isolate-specific override. Always return the control
 // to "Use map setting" on isolate changes; otherwise a prior best_read or
 // assemble click can silently override the next isolate's own map setting.
 $('#mode').value='map';
 $('#runBtnSidebar').disabled=false;$('#emptyHint').classList.add('hidden');$('#samplePanel').classList.remove('hidden');$('#sampleHeading').textContent=id;$('#sampleSub').textContent=`${s.read_count} read(s) · map default: ${s.processing_mode}`;
 $('#readRows').innerHTML=s.reads.map((r,i)=>`<tr data-file="${esc(r.file)}"><td title="${esc(r.file)}">${esc(r.file.split('/').pop())}<div class="small">${esc(r.file)}</div></td><td><input class="primer" value="${esc(r.primer)}" placeholder="e.g. 27F"></td><td><select class="direction"><option value="">Choose…</option><option ${r.direction==='forward'?'selected':''} value="forward">forward</option><option ${r.direction==='reverse'?'selected':''} value="reverse">reverse</option></select></td><td>${esc(r.processing_mode||'map default')}</td></tr>`).join('');renderSamples()
}
function readSettings(){let o={};SETTING_INFO.forEach(([k])=>o[k]=Number($(`#s_${k}`).value));return o}
const SEQUENCE_PROFILES={
 '16S':{read:300,final:800,label:'16S rRNA'},
 'ITS':{read:100,final:150,label:'ITS'},
 '18S':{read:200,final:300,label:'18S rRNA'},
 '28S_LSU':{read:200,final:300,label:'28S / LSU rRNA'},
 'CUSTOM':{read:null,final:null,label:'Other / custom target'}
};
function applySequenceProfile(){
 const key=$('#sequenceType').value||'16S',p=SEQUENCE_PROFILES[key];
 if(p.read!==null){$('#s_min_read_length').value=p.read;$('#s_min_final_length').value=p.final;}
 $('#sequenceTypeHelp').textContent=key==='CUSTOM'
  ?'No thresholds are changed. Set minimum read and final lengths for the assay below.'
  :`${p.label} starting thresholds: minimum read ${p.read} bp; minimum final sequence ${p.final} bp. These are editable starting points, not biological limits.`;
}
async function applyLoadedMap(data){
 if(data.cancelled)return;
 studioState.samples=data.samples||[];studioState.map_path=data.path||'';studioState.input_source=data.input_kind||studioState.input_source;sampleRows=data.samples||[];
 if(data.path)$('#mapPath').value=data.path;else if(data.quick_map)$('#mapPath').value='';
 if(data.taxonomy_loaded===false||['standalone_sample_map','normalised_read_map','raw_chromatograms'].includes(data.input_kind)){
  studioState.onboarding_taxonomy={loaded:false,message:'No pre-screened taxonomy is loaded (NOT_SCREENED).'};
  cohortReview=cohortReview||null;
 }
 $('#mapStatus').textContent=data.status_message||`Loaded ${studioState.samples.length} isolate(s), ${data.read_count||0} read(s).`;
 renderOnboardingTaxonomyState(studioState.onboarding_taxonomy);
 activeId='';lastInspection=null;$('#samplePanel').classList.add('hidden');$('#emptyHint').classList.remove('hidden');$('#results').classList.add('hidden');renderSamples();
 if(studioState.samples.length)selectSample(studioState.samples[0].sequence_id);
}
async function loadMap(path){setStatus($('#mapStatus'),'Loading map and checking referenced files…');try{let data=await api('/api/load-map',{path});await applyLoadedMap(data)}catch(e){setStatus($('#mapStatus'),e.message,true)}}
async function chooseMap(){setStatus($('#mapStatus'),'Choose a sample map in the file dialog…');try{let data=await api('/api/pick-map',{});await applyLoadedMap(data);if(!data.cancelled)$('#mapStatus').textContent+= ' All read paths were resolved.'}catch(e){setStatus($('#mapStatus'),e.message,true)}}
async function chooseReads(){setStatus($('#mapStatus'),'Choose AB1/ABI files in the file dialog…');try{let data=await api('/api/pick-reads',{});data.quick_map=true;data.status_message=`Loaded ${data.samples.length} isolate(s), ${data.read_count} read(s) from selected chromatograms.`;await applyLoadedMap(data);if(!data.cancelled)$('#mapStatus').textContent+=' Check inferred IDs, primers and directions before processing.'}catch(e){setStatus($('#mapStatus'),e.message,true)}}
async function chooseInputFolder(){
 const button=$('#chooseInputFolderBtn');button.disabled=true;
 setStatus($('#mapStatus'),'Choose a folder. Studio is checking for Assistant cohort files, sample maps and AB1/ABI chromatograms…');
 try{
  const d=await api('/api/pick-input-folder',{});if(d.cancelled)return;
  studioState.onboarding_taxonomy=d.onboarding_taxonomy||{loaded:Boolean(d.taxonomy_loaded),root:'',read_count:0,taxonomy_read_count:0,hit_count:0,isolate_count:0,message:d.status_message||'No pre-screened taxonomy loaded.'};
  if(d.screen_ref){studioState.screen_ref=d.screen_ref;$('#screenRefPath').value=d.screen_ref;}
  if(d.screen_taxa!==undefined){studioState.screen_taxa=d.screen_taxa||'';$('#screenTaxaPath').value=studioState.screen_taxa;}
  await applyLoadedMap(d);
  renderOnboardingTaxonomyState(studioState.onboarding_taxonomy);
  if(d.taxonomy_loaded){
   updateTaxonomyScreenStatus('Pre-screened Assistant taxonomy loaded. VSEARCH remains optional.');
  }else{
   const refAvailable=Boolean(studioState.screen_ref);
   updateTaxonomyScreenStatus(refAvailable?'No cached taxonomy (NOT_SCREENED). Fresh VSEARCH is available but opt-in.':'No cached taxonomy (NOT_SCREENED). Choose a reference FASTA to enable optional VSEARCH.');
   $('#sampleTaxonomyStatus').textContent='NOT_SCREENED';
  }
 }catch(e){setStatus($('#mapStatus'),e.message,true)}
 finally{button.disabled=false}
}
async function chooseOutput(){try{let d=await api('/api/pick-output',{});if(d.cancelled)return;studioState.out=d.out;$('#outputPath').textContent=d.out;$('#outputBadge').textContent='Output: '+d.out}catch(e){setStatus($('#batchStatus'),e.message,true)}}
function updateTaxonomyScreenStatus(message){$('#taxScreenStatus').textContent=message;$('#taxScreenConfigSource').textContent=studioState.taxonomy_config_path?`Auto-detected from ${studioState.taxonomy_config_path}`:''}
function renderOnboardingTaxonomyState(data){
 const loaded=Boolean(data&&data.loaded);
 $('#onboardingTaxonomyStatus').textContent=loaded?'Pre-screened loaded':'NOT_SCREENED';
 $('#onboardingTaxonomyDetails').textContent=loaded
  ? `${data.isolate_count} isolates · ${data.read_count} mapped reads · ${data.taxonomy_read_count} taxonomy assignments · ${data.hit_count} candidate hits · ${data.partner_count||data.map_paths?.length||0} paired dataset(s)${data.partner_names?.length?` · ${data.partner_names.join(', ')}`:''}`
  : 'No pre-screened taxonomy is loaded. This is valid for standalone Sanger work. Select a reference FASTA and opt into VSEARCH when fresh taxonomy is needed.';
}
async function applyTaxonomySettings(){
 try{let d=await api('/api/set-taxonomy-screen',{screen_ref:$('#screenRefPath').value,screen_taxa:$('#screenTaxaPath').value});studioState.screen_ref=d.screen_ref;studioState.screen_taxa=d.screen_taxa;studioState.taxonomy_screen_enabled=d.enabled;updateTaxonomyScreenStatus(d.enabled?'Reference saved. Fresh VSEARCH runs only when the re-screen checkbox is enabled; otherwise cached onboarding taxonomy is used when available.':'No reference configured; VSEARCH re-screening unavailable.');}
 catch(e){setStatus($('#taxScreenStatus'),e.message,true)}
}
async function chooseTaxonomyRef(){try{let d=await api('/api/pick-taxonomy-reference',{});if(d.cancelled)return;$('#screenRefPath').value=d.screen_ref;$('#screenTaxaPath').value=d.screen_taxa;studioState.screen_ref=d.screen_ref;studioState.screen_taxa=d.screen_taxa;studioState.taxonomy_screen_enabled=d.enabled;updateTaxonomyScreenStatus('Reference selected for optional re-screening.')}catch(e){setStatus($('#taxScreenStatus'),e.message,true)}}
async function chooseTaxonomyTable(){try{let d=await api('/api/pick-taxonomy-table',{});if(d.cancelled)return;$('#screenRefPath').value=d.screen_ref;$('#screenTaxaPath').value=d.screen_taxa;studioState.screen_ref=d.screen_ref;studioState.screen_taxa=d.screen_taxa;studioState.taxonomy_screen_enabled=d.enabled;updateTaxonomyScreenStatus('Taxonomy mapping selected for optional re-screening.')}catch(e){setStatus($('#taxScreenStatus'),e.message,true)}}
async function shutdownStudio(){
 const button=$('#shutdownBtn');
 if(inspectionRunning||batchRunning){setStatus($('#mapStatus'),'Wait for the current inspection or batch processing run to finish before shutting down Studio.',true);return}
 if(!window.confirm('Shut down the local Sanger Studio server? Saved inspections, outputs and review decisions will remain on disk.'))return;
 button.disabled=true;button.textContent='Shutting down…';
 try{
  const result=await api('/api/shutdown',{});
  $('#shutdownNotice').textContent=result.message||'Sanger Studio has stopped. Your saved outputs and review decisions remain on disk. You can close this browser tab.';
  $('#shutdownNotice').classList.add('show');
  $('#mapStatus').textContent='Sanger Studio server stopped. You can close this browser tab.';
  $('#outputBadge').textContent='Studio stopped';
  document.querySelectorAll('button').forEach(b=>{if(b!==button)b.disabled=true});
  button.textContent='Studio stopped';
 }catch(e){button.disabled=false;button.textContent='Shut down Studio';setStatus($('#mapStatus'),e.message,true)}
}
async function attachReview(){setStatus($('#cohortReviewStatus'),'Choose a cohort marker review template…');try{let d=await api('/api/pick-marker-review',{});if(d.cancelled)return;renderCohortReview(d);$('#cohortReviewStatus').textContent=`Attached ${d.source} · ${d.total} marker(s).`; }catch(e){setStatus($('#cohortReviewStatus'),e.message,true)}}
function renderCohortReview(data){cohortReview=data;$('#cohortReviewControls').classList.remove('hidden');$('#cohortReviewPill').className='pill neutral';$('#cohortReviewPill').textContent=`${data.decided}/${data.total} decided`;$('#cohortReviewStatus').textContent=`Source: ${data.source}`;$('#cohortReviewer').value=$('#cohortReviewer').value||((data.rows.find(r=>r.reviewer)?.reviewer)||'');renderCohortReviewRows();}
function renderReviewSequenceEvidence(r){
 const rawLength=String(r.output_length??'').trim();
 const outputLength=Number(rawLength);
 const hasFinal=rawLength!==''&&Number.isFinite(outputLength)&&outputLength>0;
 const nValue=String(r.n_percent??'').trim();
 const outputType=String(r.output_type||'').toUpperCase();
 const outputLabel=outputType==='INDIVIDUAL_READ'?'Selected individual read':'Final consensus';
 const finalLine=hasFinal
  ? `<div><b>${esc(outputLabel)}:</b> ${esc(outputLength)} bp · N: ${nValue?`${esc(nValue)}%`:'not recorded'}</div>`
  : `<div><b>Final consensus:</b> none (${rawLength||'0'} bp) · N: not applicable</div>`;
 const candidateLength=Number(r.candidate_length??0);
 const candidateSource=String(r.candidate_source??'').trim();
 const candidateN=String(r.candidate_n_percent??'').trim();
 const candidateRead=String(r.candidate_read_id??'').trim();
 let candidateLine='';
 if(candidateLength>0&&candidateSource==='BEST_RECOVERABLE_READ'){
  candidateLine=`<div><b>Fallback read candidate:</b> ${esc(candidateLength)} bp · N: ${candidateN?`${esc(candidateN)}%`:'not recorded'}${candidateRead?`<div class="small">Read: ${esc(candidateRead)} · not an assembled consensus</div>`:'<div class="small">Not an assembled consensus</div>'}</div>`;
 }else if(candidateLength>0&&candidateSource==='FAILED_FINAL_CONSENSUS'){
  candidateLine=`<div><b>Review candidate:</b> failed final consensus (${esc(candidateLength)} bp · N: ${candidateN?`${esc(candidateN)}%`:'not recorded'})</div>`;
 }
 const readEvidence=String(r.read_level_evidence??'').trim();
 const readLine=readEvidence?`<div class="reviewReadEvidence"><b>Read-level QC:</b> ${esc(readEvidence)}</div>`:'';
 return finalLine+candidateLine+readLine;
}
function renderFallbackReadChoice(r){
 const candidates=Array.isArray(r.fallback_read_options_list)?r.fallback_read_options_list:[];
 if(candidates.length<2)return '<span class="taxon-muted">Not required</span>';
 const recommended=String(r.recommended_read_id||candidates[0].read_id||'');
 const selected=String(r.selected_read_id||recommended);
 const source=String(r.read_selection_source||'').toUpperCase();
 const options=candidates.map((candidate,i)=>{
  const id=String(candidate.read_id||'');
  const label=`${i===0?'Recommended · ':''}${id} · ${candidate.length||0} bp · Q${candidate.mean_quality??'NA'} · N ${candidate.n_percent??'NA'}%`;
  return `<option value="${esc(id)}" ${id===selected?'selected':''}>${esc(label)}</option>`;
 }).join('');
 return `<select class="fallbackReadChoice" aria-label="Fallback read for ${esc(r.sequence_id)}">${options}</select><div class="small">${source==='MANUAL'?'Manual reviewer selection':'Automatic recommendation: longest usable read'}</div>`;
}
function renderCohortReviewRows(){
 if(!cohortReview)return;
 let q=$('#cohortReviewSearch').value.toLowerCase(),filter=$('#cohortReviewFilter').value;
 let rows=cohortReview.rows.filter(r=>(filter==='all'||String(r.qc_class).toUpperCase()===filter)&&[r.sequence_id,r.partner_id,r.dataset,r.qc_reasons,r.taxonomy_concordance,r.read_level_evidence,r.candidate_source,r.candidate_read_id,(r.fallback_read_options_list||[]).map(x=>x.read_id).join(' ')].join(' ').toLowerCase().includes(q));
 $('#cohortReviewRows').innerHTML=rows.map(r=>{
  let options=['<option value="">— Not reviewed —</option>',...(r.allowed_decisions_list||[]).map(d=>`<option value="${esc(d)}" ${String(r.decision||'').toUpperCase()===d?'selected':''}>${esc(d)}</option>`)].join('');
  let cls=String(r.qc_class).toUpperCase()==='FAIL_QC'?'qcfail':String(r.qc_class).toUpperCase()==='PASS_WITH_WARNINGS'?'warning':'';
  const decision=String(r.decision||'').toUpperCase();
  const overrideHint=decision==='APPROVE_FAIL_QC'?'<div class="small review-override-hint">Manual override recorded · automated QC remains FAIL_QC</div>':'';
  return `<tr class="${cls}" data-seq="${esc(r.sequence_id)}"><td><button class="quiet reviewInspect" data-review-inspect="${esc(r.sequence_id)}">${esc(r.sequence_id)}</button><div class="small">${esc(r.marker_source_type||'')}</div></td><td>${esc(r.partner_id||'')}<div class="small">${esc(r.dataset||'')}</div></td><td><span class="pill ${pillClass(r.qc_class)}">${esc(r.qc_class)}</span><div class="small">${esc(r.review_type||'')}</div><div class="small">${esc(r.processing_mode||'')} · ${esc(r.assembly_status||'')}</div></td><td>${esc(r.qc_reasons||'—')}</td><td>${renderReviewSequenceEvidence(r)}</td><td>${renderFallbackReadChoice(r)}</td><td>${esc(r.taxonomy_concordance||'—')}<div class="small">${esc(r.read_taxonomy_assignments||'')}</div></td><td><select class="cohortDecision">${options}</select><div class="small">${(r.allowed_decisions_list||[]).join(' · ')}</div>${overrideHint}</td><td><input class="cohortNotes" type="text" value="${esc(r.notes||'')}" placeholder="Reason / evidence"></td></tr>`;
 }).join('')||'<tr><td colspan="9">No records match this filter.</td></tr>';
 let counts={};cohortReview.rows.forEach(r=>{let key=String(r.qc_class||'UNKNOWN');counts[key]=(counts[key]||0)+1});
 let decided=cohortReview.rows.filter(r=>String(r.decision||'').trim()).length;
 let manualSelections=cohortReview.rows.filter(r=>String(r.read_selection_source||'').toUpperCase()==='MANUAL'&&String(r.selected_read_id||'').trim()).length;
 $('#cohortReviewCounts').textContent=`Showing ${rows.length} of ${cohortReview.total} · ${decided} decisions · ${manualSelections} manual read selection(s) · `+Object.entries(counts).map(([k,v])=>`${k}: ${v}`).join(' · ');
 $('#cohortReviewPill').textContent=`${decided}/${cohortReview.total} decided`;
}

async function saveCohortReview(){
 if(!cohortReview)return;
 let rows=[...document.querySelectorAll('#cohortReviewRows tr[data-seq]')].map(tr=>({
  sequence_id:tr.dataset.seq,
  decision:tr.querySelector('.cohortDecision').value,
  notes:tr.querySelector('.cohortNotes').value,
  selected_read_id:tr.querySelector('.fallbackReadChoice')?.value||''
 }));
 /* Include filtered-out rows unchanged so selections and decisions are preserved. */
 let visible=new Map(rows.map(r=>[r.sequence_id,r]));
 let all=cohortReview.rows.map(r=>visible.get(r.sequence_id)||({
  sequence_id:r.sequence_id,decision:r.decision||'',notes:r.notes||'',
  selected_read_id:r.selected_read_id||''
 }));
 const reviewer=$('#cohortReviewer').value.trim();
 const metadata=new Map(cohortReview.rows.map(r=>[r.sequence_id,r]));
 const problems=[];
 for(const row of all){
  const meta=metadata.get(row.sequence_id)||{};
  const decision=String(row.decision||'').trim().toUpperCase();
  const rowReviewer=reviewer||String(meta.reviewer||'').trim();
  if(decision&&!rowReviewer){
   problems.push(`${row.sequence_id}: enter a reviewer name or initials.`);
  }
  if(decision.startsWith('APPROVE')&&!String(row.notes||'').trim()){
   problems.push(`${row.sequence_id}: an approval needs explanatory notes.`);
  }
  if(decision&&!(meta.allowed_decisions_list||[]).includes(decision)){
   problems.push(`${row.sequence_id}: ${decision} is not an allowed decision for this QC class.`);
  }
  const selected=String(row.selected_read_id||'').trim();
  const options=meta.fallback_read_options_list||[];
  const optionIds=options.map(x=>String(x.read_id||'').trim());
  if(selected&& !optionIds.includes(selected)){
   problems.push(`${row.sequence_id}: selected fallback read is not in the QC-passing options.`);
  }
  if(selected&&selected!==String(meta.recommended_read_id||'').trim()&&!rowReviewer){
   problems.push(`${row.sequence_id}: a manual fallback-read choice needs reviewer identity.`);
  }
 }
 if(problems.length){
  const visibleProblems=problems.slice(0,8);
  const remainder=problems.length-visibleProblems.length;
  const lineBreak=String.fromCharCode(10);
  const message=`Cannot save yet (${problems.length} issue(s)):${lineBreak}${visibleProblems.join(lineBreak)}${remainder?`${lineBreak}…and ${remainder} more. Filter the table and correct the remaining records.`:''}`;
  setStatus($('#cohortReviewSaveStatus'),message,true);
  $('#cohortReviewSaveStatus').scrollIntoView({behavior:'smooth',block:'nearest'});
  return;
 }
 let button=$('#saveCohortReviewBtn');button.disabled=true;
 setStatus($('#cohortReviewSaveStatus'),`Saving ${all.length} marker review row(s)…`);
 try{
  let d=await api('/api/save-marker-review',{reviewer,rows:all});
  cohortReview=d;renderCohortReview(d);
  $('#cohortReviewSaveStatus').innerHTML=`Saved ${d.decided} decision(s) and manual read selections to <b>${esc(d.saved)}</b>. <a href="${esc(d.download_url)}" download>Download completed review TSV</a>. Supply this file to Assistant PREPARE with --marker-review to carry the choices forward.`;
 }catch(e){
  setStatus($('#cohortReviewSaveStatus'),`Save failed: ${e.message}`,true);
  $('#cohortReviewSaveStatus').scrollIntoView({behavior:'smooth',block:'nearest'});
 }finally{button.disabled=false}
}

function setRunFeedback(message,{error=false,loading=false}={}){
 ['#runStatus','#sidebarRunStatus'].forEach(sel=>{
  const node=$(sel);if(!node)return;
  node.className=error?'error':'statusline small';
  node.innerHTML=loading?`<span class="spinner"></span> ${esc(message)}`:esc(message);
 });
}
async function inspect(){
 if(!activeId||inspectionRunning)return;
 const inspectedId=activeId;
 const requestedMode=$('#mode').value||'map';
 const sample=(studioState.samples||[]).find(x=>x.sequence_id===inspectedId);
 const mapMode=sample?.processing_mode||'assemble';
 const effectiveMode=requestedMode==='map'?mapMode:requestedMode;
 let rows=[...document.querySelectorAll('#readRows tr')].map(r=>({file:r.dataset.file,primer:r.querySelector('.primer').value,direction:r.querySelector('.direction').value}));
 if(rows.some(r=>!r.primer||!r.direction)){setRunFeedback('Set a primer and direction for every read before processing.',{error:true});return}
 const rerun=$('#rerunTaxonomyCheck').checked;
 const reference=$('#screenRefPath').value||studioState.screen_ref||'configured reference';
 const referenceName=String(reference).split(/[\\/]/).pop();
 const modeLabel=effectiveMode==='best_read'?'single best read (assembly skipped)':'assemble primer reads';
 const body={sequence_id:inspectedId,reads:rows,processing_mode:requestedMode,settings:readSettings(),sequence_type:$('#sequenceType').value,custom_primer_sequences:$('#customPrimerSequences').value,rerun_taxonomy:rerun};
 inspectionRunning=true;$('#shutdownBtn').disabled=true;
 const runButtons=[$('#runBtn'),$('#runBtnSidebar')].filter(Boolean);
 const originalButtonLabels=runButtons.map(b=>b.textContent);
 runButtons.forEach(b=>{b.disabled=true;b.innerHTML=rerun?'<span class="spinner"></span> Inspect + VSEARCH running…':'<span class="spinner"></span> Inspecting isolate…'});
 const runningMessage=rerun
  ? `Inspecting ${inspectedId} using ${modeLabel} and running a FRESH VSEARCH taxonomy screen against ${referenceName} on the current trimmed reads. VSEARCH is running as part of this inspection and may take several minutes…`
  : `Inspecting ${inspectedId} using ${modeLabel}; taxonomy will use pre-screened onboarding results where available…`;
 setRunFeedback(runningMessage,{loading:true});
 if(rerun)updateTaxonomyScreenStatus(`RUNNING: fresh VSEARCH for ${inspectedId} against ${referenceName} using the current trimmed reads…`);
 $('#results').classList.add('hidden');
 try{
  let r=await api('/api/inspect',body);lastInspection=r;
  // The backend stores an explicit per-isolate override in this Studio session.
  // Keep the client-side map summary in sync so returning to this isolate and
  // choosing "Use map setting" reproduces the last resolved mode for this row.
  const resolvedSample=(studioState.samples||[]).find(x=>x.sequence_id===inspectedId);
  if(resolvedSample&&r.processing_mode){
   resolvedSample.processing_mode=r.processing_mode;
   (resolvedSample.reads||[]).forEach(read=>{read.processing_mode=r.processing_mode});
  }
  renderResult(r);
  setRunFeedback(rerun?(r.taxonomy_screen_message||'Fresh VSEARCH taxonomy screening completed as part of this inspection.'):'Inspection complete.');
  if(rerun)updateTaxonomyScreenStatus(r.taxonomy_screen_message||'Fresh VSEARCH completed.');
 }catch(e){
  setRunFeedback(rerun?`Inspection/VSEARCH failed: ${e.message}`:e.message,{error:true});
  if(rerun)updateTaxonomyScreenStatus(`FAILED: ${e.message}`);
 }finally{
  inspectionRunning=false;$('#shutdownBtn').disabled=false;
  runButtons.forEach((b,i)=>{b.disabled=!activeId;b.textContent=originalButtonLabels[i]||'Inspect and assemble this isolate'});
 }
}
async function processAll(){
 let rows=[...document.querySelectorAll('#readRows tr')].map(r=>({file:r.dataset.file,primer:r.querySelector('.primer').value,direction:r.querySelector('.direction').value}));
 let rerun=$('#rerunTaxonomyCheck').checked,reference=$('#screenRefPath').value||studioState.screen_ref||'configured reference';
 let body={sequence_id:activeId,reads:rows,processing_mode:$('#mode')?.value||'map',settings:readSettings(),sequence_type:$('#sequenceType').value,custom_primer_sequences:$('#customPrimerSequences').value,rerun_taxonomy:rerun};
 let button=$('#processAllBtn');button.disabled=true;button.innerHTML='<span class="spinner"></span> Processing all mapped isolates…';
 batchRunning=true;$('#shutdownBtn').disabled=true;
 setStatus($('#batchStatus'),rerun?`Processing the batch and running fresh VSEARCH taxonomy screening against ${String(reference).split(/[\\/]/).pop()} using current trimmed reads. This may take several minutes…`:'Running the shared Paper Trail/Merge Meeting pipeline for the full map. This may take a while.');
 if(rerun)updateTaxonomyScreenStatus(`RUNNING: fresh VSEARCH batch screening against ${String(reference).split(/[\\/]/).pop()}…`);
 try{
  let r=await api('/api/process-all',body);let counts=Object.entries(r.qc_counts||{}).map(([k,v])=>`${k}: ${v}`).join(' · ');
  let tax=r.taxonomy_mode==='fresh_vsearch'?` · ${r.taxonomy_screen_message||'Fresh VSEARCH re-screen completed'}`:(r.taxonomy_mode==='assistant_onboarding_cache'?' · reused onboarding taxonomy':'');
  let freshness=r.taxonomy_freshness==='STALE_TRIMMED_SEQUENCE_CHANGED'?' · WARNING: cached taxonomy is stale for changed trimmed reads':(r.taxonomy_freshness==='MATCHES_CURRENT_TRIMMED_READS'?' · cached taxonomy matches current trimmed reads':'');
  $('#batchStatus').innerHTML=`Processed ${r.sample_count} isolate(s). ${esc(counts)}${esc(tax)}${esc(freshness)}. <a href="${esc(r.bundle_url)}" download>Download complete results bundle</a>`;
  if(rerun)updateTaxonomyScreenStatus(r.taxonomy_screen_message||'Fresh VSEARCH batch screening completed.');
  let link=$('#batchStatus a');link.click();
 }catch(e){
  setStatus($('#batchStatus'),rerun?`Batch/VSEARCH failed: ${e.message}`:e.message,true);
  if(rerun)updateTaxonomyScreenStatus(`FAILED: ${e.message}`);
 }
 finally{batchRunning=false;$('#shutdownBtn').disabled=false;button.disabled=false;button.textContent='Process and download all'}
}

function renderTaxonomyResults(r){
 const host=$('#taxonomyResults');
 if(!r.taxonomy_screen_requested&&r.taxonomy_freshness!=='FRESH_VSEARCH_RERUN'){
  host.innerHTML='<div class="taxonomy-results"><h3>Primer-read taxonomy evidence</h3><p class="taxon-muted">No pre-screened taxonomy assignments were found for these read IDs. This is normal for standalone input: the status is NOT_SCREENED. Configure a reference FASTA and enable optional VSEARCH to generate fresh assignments.</p></div>';
  return;
 }
 const reads=r.taxonomy_reads||[],hits=r.taxonomy_hits||[];
 const concordance=r.taxonomy_concordance||'NOT_REPORTED';
 let html=`<div class="taxonomy-results"><h3>Primer-read taxonomy evidence</h3><p><b>Source:</b> ${esc(r.taxonomy_source||'Taxonomy screen')} · <b>Inter-read comparison:</b> ${esc(concordance)} <span class="taxon-muted">${esc(r.read_taxonomy_assignments||'')}</span></p>`;
 if(r.taxonomy_screen_message)html+=`<p class="taxon-muted"><b>${esc(r.taxonomy_freshness==='FRESH_VSEARCH_RERUN'?'Fresh VSEARCH run: ':'Taxonomy screen: ')}</b>${esc(r.taxonomy_screen_message)}</p>`;
 if(r.taxonomy_staleness_note)html+=`<p class="taxon-muted">${esc(r.taxonomy_staleness_note)}</p>`;
 if(r.taxonomy_comparison_reason)html+=`<p class="taxon-muted">${esc(r.taxonomy_comparison_reason)}</p>`;
 if(!reads.length){
  html+='<p class="taxon-muted">The taxonomy screen ran, but no per-read assignments were produced. Check the taxonomy-screen output files and the run manifest.</p>';
 } else {
  html+='<h4>Per-read taxonomy assignments</h4><div class="readscroll"><table class="reads qualityreads taxonomy-summary"><thead><tr><th>Read</th><th>Best reference hit</th><th>Identity</th><th>Query coverage</th><th>Target coverage</th><th>Taxonomy</th><th>Assigned rank</th><th>Qualified hits</th><th>Near-best hits</th><th>Classification status</th><th>Support</th></tr></thead><tbody>';
  html+=reads.map(x=>`<tr><td>${esc(x.ID)}</td><td>${esc(x.BestHit||'—')}</td><td>${esc(x.Identity||'—')}${x.Identity?'%':''}</td><td>${esc(x.QueryCoverage||'—')}${x.QueryCoverage?'%':''}</td><td>${esc(x.TargetCoverage||'—')}${x.TargetCoverage?'%':''}</td><td>${esc(x.Taxon||'No taxonomy assigned')}</td><td>${esc(x.AssignedRank||'NA')}</td><td>${esc(x.QualifiedHitCount??'0')}</td><td>${esc(x.NearBestHitCount??'0')}</td><td>${esc(x.ReadClassificationStatus||'UNKNOWN')}</td><td>${esc(x.TaxonomySupport||'—')}</td></tr>`).join('');
  html+='</tbody></table></div>';
  html+='<h4>Candidate reference hits</h4><p class="taxon-muted">Near-best hits are within the configured identity window of the best qualifying hit. Qualified indicates that identity and query/target coverage filters were met.</p>';
  html+='<div class="readscroll"><table class="reads qualityreads taxonomy-hits"><thead><tr><th>Read</th><th>Reference ID</th><th>Taxonomy</th><th>Identity</th><th>Query coverage</th><th>Target coverage</th><th>Alignment (bp)</th><th>Qualified</th><th>Near-best</th><th>Delta from best</th></tr></thead><tbody>';
  html+=hits.map(x=>`<tr><td>${esc(x.ReadID)}</td><td>${esc(x.TargetID)}</td><td>${esc(x.Taxonomy||'Unmapped reference')}</td><td>${esc(x.Identity||'—')}%</td><td>${esc(x.QueryCoverage||'—')}%</td><td>${esc(x.TargetCoverage||'—')}%</td><td>${esc(x.AlignmentLength||'—')}</td><td>${String(x.Qualified||'').toLowerCase()==='yes'?'<span class="taxon-qualified">Yes</span>':'No'}</td><td>${String(x.NearBest||'').toLowerCase()==='yes'?'<span class="taxon-nearbest">Near-best</span>':'No'}</td><td>${esc(x.IdentityDelta||'—')}</td></tr>`).join('')||'<tr><td colspan="10">No qualifying reference hits were recovered for these reads.</td></tr>';
  html+='</tbody></table></div>';
 }
 html+='</div>';host.innerHTML=html;
}
function renderResult(r){
 $('#results').classList.remove('hidden');$('#resultHeading').textContent=r.sequence_id;$('#resultPill').className='pill '+pillClass(r.qc_class);$('#resultPill').textContent=`${r.qc_class} · ${r.recommendation}`;$('#jobDir').textContent=r.job_dir;
 let ms=[['Result',r.status],['Processing mode',r.processing_mode||'map default'],['Consensus / selected read',r.output_length?`${r.output_length} bp`:'none'],['Reads used',`${r.used_reads}/${r.read_count}`],['Overlap',r.overlap_length?`${r.overlap_length} bp`:'none'],['Overlap identity',r.overlap_identity||'NA']];$('#metrics').innerHTML=ms.map(x=>`<div class="metric"><label>${esc(x[0])}</label><b>${esc(x[1])}</b></div>`).join('');
 let selectionLine=r.selected_read_id?`<br><b>Selected individual read:</b> ${esc(r.selected_read_id)} (${esc(r.read_selection_source||'selection')})${r.fallback_reason?` · ${esc(r.fallback_reason)}`:''}`:'';
 let fallbackLine=r.fallback_read_options?.length>1?`<br><b>Fallback candidates:</b> ${r.fallback_read_options.map(x=>`${esc(x.read_id)} (${esc(x.length)} bp, Q${esc(x.mean_quality)}, N ${esc(x.n_percent)}%)`).join(' · ')}<br><span class="taxon-muted">Choose the read to carry forward in the Fallback read selection column of the Cohort marker review table, then save the review TSV.</span>`:'';
 $('#why').innerHTML=`<b>Processing mode:</b> ${esc(r.processing_mode||'map default')} · <b>Assembly outcome:</b> ${esc(r.assembly_status||r.status)}<br><b>QC reason:</b> ${esc(r.reasons||'No QC warning or failure reason.')}<br><b>Suggested action:</b> ${esc(r.suggested_action||'Review the evidence and retain the report with the submission.')}${selectionLine}${fallbackLine}${r.taxonomy_concordance?`<br><b>Primer-read taxonomy agreement:</b> ${esc(r.taxonomy_concordance)} · ${esc(r.read_taxonomy_assignments)}`:''}`;
 if(cohortReview){
  const reviewRow=cohortReview.rows.find(x=>x.sequence_id===r.sequence_id);
  if(reviewRow&&Array.isArray(r.fallback_read_options)&&r.fallback_read_options.length>1){
   reviewRow.fallback_read_options_list=r.fallback_read_options;
   reviewRow.fallback_read_options=JSON.stringify(r.fallback_read_options);
   reviewRow.recommended_read_id=r.recommended_read_id||r.fallback_read_options[0].read_id;
   reviewRow.selected_read_id=r.selected_read_id||reviewRow.recommended_read_id;
   reviewRow.read_selection_source=r.read_selection_source||'AUTO_LONGEST_USABLE_READ';
   reviewRow.fallback_reason=r.fallback_reason||'';
   renderCohortReviewRows();
  }
 }

 const ev=r.consensus_evidence||[];
 const consensus=ev.map(x=>String(x.ConsensusBase||'N')).join('');
 const disagreements=ev.filter(x=>String(x.Disagreement).toLowerCase()==='yes');
 let activeEvidence=ev;
 let activeDisagreements=disagreements;
 let activeReviewPositions=[]; // candidate-specific disagreements plus ambiguous-call review sites
 const single=ev.filter(x=>String(x.EvidenceClass||'').includes('SINGLE_READ')).length;
 const candidates=(r.assembly_candidates||[]).map(x=>({...x, len:Number(x.OutputLength)||0}));
 const usable=candidates.filter(x=>x.len>0&&!String(x.Reasons||'').includes('partial_no_overlap'));
 const ranked=usable.sort((a,b)=>(a.QCClass==='PASS_HIGH_CONFIDENCE'?-1:0)-(b.QCClass==='PASS_HIGH_CONFIDENCE'?-1:0)||b.len-a.len);
 const fasta=`>${r.sequence_id}|SangerStudio_consensus|${r.qc_class}\n${(consensus.match(/.{1,70}/g)||[]).join('\n')}\n`;
 const dl=(name,content,type='text/plain')=>{const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([content],{type}));a.download=name;a.click();URL.revokeObjectURL(a.href)};
 const pct=ev.length?Math.round(100*disagreements.length/ev.length):0;
 const seqHtml=consensus?`<div class="consensus-seq" id="consensusSeq">${[...consensus].map((b,i)=>`<span class="base ${b==='N'?'ambiguous':''} ${String(ev[i]?.Disagreement).toLowerCase()==='yes'?'discordant':''}" title="Position ${i+1}; ${esc(ev[i]?.ObservedBases||'')}; ${esc(ev[i]?.EvidenceClass||'')}">${b}</span>`).join('')}</div>`:'<p class="muted">No consensus-position evidence was supplied for this result.</p>';
 $('#readResults').innerHTML=`
 <section class="review-hero">
   <div><div class="eyebrow">SEQUENCE REVIEW · ${esc(r.sequence_id)}</div><h2>Assembly overview</h2><p class="muted">Review the sequence first. Supporting evidence and technical reports are available below.</p></div>
   <div class="hero-status ${r.qc_class==='FAIL_QC'?'fail':r.qc_class==='PASS_WITH_WARNINGS'?'warn':'pass'}">${esc(r.qc_class||r.status||'UNCLASSIFIED')}</div>
 </section>
 <div class="review-metrics">
  <div><b>${Number(r.output_length)||consensus.length||'—'}</b><span>Consensus bp</span></div>
  <div><b>${esc(r.used_reads||'—')}</b><span>Reads used</span></div>
  <div><b>${disagreements.length}</b><span>Disagreement positions</span></div>
  <div><b>${single}</b><span>Single-read positions</span></div>
 </div>
 <section class="review-panel">
  <div class="panel-heading"><div><h3>Consensus sequence</h3><p>Ambiguous bases and read disagreements are highlighted. Select a base to inspect its evidence.</p></div><button class="blue" id="downloadConsensus">Download consensus FASTA</button></div>
  <div class="sequence-toolbar"><span><i class="legend-dot blue-dot"></i> Consensus</span><span><i class="legend-dot red-dot"></i> Read disagreement</span><span><i class="legend-dot amber-dot"></i> Ambiguous base</span><span class="muted">${consensus.length.toLocaleString()} positions</span></div>
  <div class="sequence-scroll">${seqHtml}</div>
  <div class="sequence-inspector" id="baseInspector">Select a nucleotide above to see its position-level evidence.</div>
  <div class="review-note">${disagreements.length} disagreement positions (${pct}% of positions with evidence); ${single} positions have single-read support. These are review indicators, not automatic rejection criteria.</div>
 </section>
 <section class="review-panel" id="readAlignmentPanel">
  <div class="panel-heading"><div><h3>Multiple-read alignment</h3><p>Reads are shown in the selected candidate's coordinates. 'Next disagreement' cycles through read-review positions: different unambiguous bases and ambiguous-vs-concrete calls. The candidate-card 'ambiguous overlap conflicts' count is a separate assembly/QC metric.</p></div><span class="small">${(r.read_alignment||[]).length} placed reads</span></div>
  <div class="alignment-controls"><label>Consensus position <input id="alignmentStart" type="number" min="1" max="${Math.max(1,consensus.length)}" value="1"></label><input id="alignmentRange" type="range" min="1" max="${Math.max(1,consensus.length)}" value="1"><span id="alignmentWindowLabel">1–${Math.min(80,consensus.length)} / ${consensus.length}</span><button class="secondary" id="nextDisagreement">Next disagreement</button></div>
  <div class="alignment-legend"><span class="aln-match">A/C/G/T</span> mapped base <span class="aln-gap">–</span> no mapped base <span class="aln-diff">Mismatch</span> differs from consensus <span class="aln-ambig">N</span> ambiguous call</div>
  <div class="alignment-scroll" id="readAlignment"></div>
  <p class="help" id="alignmentStatus">${(r.read_alignment||[]).length ? 'Select a position in the alignment to inspect read-level evidence.' : 'No read placements with corresponding trimmed sequences were available. Re-run inspection with the updated version to generate the alignment.'}</p>
 </section>
 <section class="review-panel">
  <div class="panel-heading"><div><h3>Select the assembly to retain</h3><p>Choose one sequence for this isolate. QC-failed reads are included only in clearly flagged review candidates; selecting one is an explicit override, not a change to automated QC.</p></div><span class="small">${ranked.length} sequence-bearing candidates</span></div>
  <div class="candidate-grid">${ranked.map((x,i)=>`<article class="candidate-card ${x.IncludesFailedRead==='yes'?'candidate-override':''}"><div class="candidate-top"><b>${esc(x.CandidateID)}</b><span class="candidate-qc ${String(x.QCClass).includes('PASS_HIGH')?'pass':String(x.QCClass).includes('PASS')?'warn':'fail'}">${esc(x.QCClass)}</span></div><div class="candidate-length">${x.len.toLocaleString()} <small>bp</small></div><p>${esc(x.ReadIDs||'Read set not recorded')}</p><div class="candidate-meta">Overlap: ${esc(x.OverlapLength||'—')} bp · Identity: ${Number.isFinite(Number(x.OverlapIdentity))?(Number(x.OverlapIdentity)*100).toFixed(2)+'%':'—'} · Ambiguous overlap conflicts: ${esc(x.AmbiguousConflicts||x.Conflicts||'0')}</div>${x.IncludesFailedRead==='yes'?`<div class="override-warning">Includes read(s) failing read-level QC: ${esc(x.ReadQC||'see read QC table')}</div>`:''}<small>${esc(x.Reasons||'No additional QC warnings reported')}</small><button class="candidate-select" data-candidate="${esc(x.CandidateID)}" ${x.Sequence?'':'disabled'}>Select this sequence</button></article>`).join('')||'<p class="muted">No sequence-bearing candidates were generated.</p>'}</div>
  <div class="selection-footer"><span id="candidateSelectionStatus">No candidate selected for final extraction.</span><button class="blue" id="downloadSelectedCandidate" disabled>Download selected FASTA</button><button class="secondary" id="downloadSelectedBatch" ${Object.keys(selectedAssemblySequences).length?'':'disabled'}>Download selected isolates FASTA (${Object.keys(selectedAssemblySequences).length})</button></div>
  <p class="help">The selected candidate is exported using its exact stored sequence from assembly_candidates.fasta. Selecting it does not alter the original automated result or QC classification.</p>
 </section>
 <details class="diagnostics"><summary>Read quality and trimming evidence (${r.reads.length} reads)</summary><div class="readscroll"><table class="reads qualityreads"><thead><tr><th>Read</th><th>Primer / direction</th><th>QC</th><th>Raw bp</th><th>Trimmed bp</th><th>Trim coordinates</th><th>Mean Q</th><th>Expected errors</th><th>Masked N%</th><th>Trace / peak evidence</th><th>Reason</th></tr></thead><tbody>${r.reads.map(x=>`<tr><td>${esc(x.ReadID)}</td><td>${esc(x.Primer)} / ${esc(x.Direction)}</td><td>${esc(x.QCClass)}</td><td>${esc(x.RawLength)}</td><td>${esc(x.TrimmedLength)}</td><td>${esc(x.TrimStart)}–${esc(x.TrimEnd)}</td><td>${esc(x.MeanTrimmedQuality)}</td><td>${esc(x.TrimmedExpectedErrors)}</td><td>${esc(x.MaskedNPercent)}%</td><td>Trace ${esc(x.TraceAvailable)} · PCON ${esc(x.QualityAvailable)} · mixed peaks ${esc(x.MixedPeakCount)}</td><td>${esc(x.Reasons||'—')}</td></tr>`).join('')}</tbody></table></div></details>
 <details class="diagnostics"><summary>Read overlaps and assembly diagnostics (${r.overlap_candidates?.length||0} pairs)</summary><div class="readscroll"><table class="reads qualityreads"><thead><tr><th>Read A</th><th>Read B</th><th>Overlap</th><th>Identity</th><th>Outcome</th><th>Reason</th></tr></thead><tbody>${(r.overlap_candidates||[]).map(x=>`<tr><td>${esc(x.ReadA)}</td><td>${esc(x.ReadB)}</td><td>${esc(x.OverlapLength)} bp</td><td>${(Number(x.Identity)*100).toFixed(2)}%</td><td>${String(x.Accepted).toLowerCase()==='yes'?'Accepted':'Rejected'}</td><td>${esc(x.Reason)}</td></tr>`).join('')}</tbody></table></div></details>
 <details class="diagnostics"><summary>Full per-position evidence (${ev.length.toLocaleString()} positions)</summary><div class="readscroll position-table"><table class="reads qualityreads"><thead><tr><th>Position</th><th>Consensus</th><th>Supporting reads</th><th>Observed bases</th><th>Disagreement</th><th>Evidence class</th></tr></thead><tbody>${ev.map(x=>`<tr><td>${esc(x.ConsensusPosition)}</td><td><b>${esc(x.ConsensusBase)}</b></td><td>${esc(x.SupportingReads)}</td><td>${esc(x.ObservedBases)}</td><td>${esc(x.Disagreement)}</td><td>${esc(x.EvidenceClass)}</td></tr>`).join('')}</tbody></table></div></details>`;
 $('#downloadConsensus').onclick=()=>{const seq=activeConsensus||consensus;const label=selectedCandidate?selectedCandidate.CandidateID:'consensus';dl(`${r.sequence_id}_${label}.fasta`,`>${r.sequence_id}|${label}|${selectedCandidate?.QCClass||r.qc_class}\n${(seq.match(/.{1,70}/g)||[]).join('\n')}\n`,'text/x-fasta');};
 let selectedCandidate=null;
 const candidateById=new Map(ranked.map(x=>[x.CandidateID,x]));
 const selectedStatus=$('#candidateSelectionStatus');
 const selectedDownload=$('#downloadSelectedCandidate');
 function activateCandidateAlignment(candidate){
   activeConsensus=String(candidate.Sequence||'').toUpperCase();
   alignmentReads=(r.candidate_read_alignments||{})[candidate.CandidateID]||[];
   // Recalculate evidence against the selected candidate, never reuse the
   // default assembly's position table for a different sequence. N/other
   // ambiguous calls are displayed but do not count as nucleotide agreement.
   activeEvidence=Array.from({length:activeConsensus.length},(_,i)=>{
     const calls=alignmentReads.map(x=>({read_id:x.read_id,base:String(x.aligned_sequence||'')[i]||'-'})).filter(x=>x.base!=='-');
     const concrete=calls.filter(x=>'ACGT'.includes(x.base));
     const bases=[...new Set(concrete.map(x=>x.base))];
     const disagreement=bases.length>1;
     return {ConsensusPosition:i+1,ConsensusBase:activeConsensus[i],SupportingReads:calls.map(x=>x.read_id).join(';'),ObservedBases:calls.map(x=>x.base).join(';'),Disagreement:disagreement?'yes':'no',EvidenceClass:disagreement?'MULTI_READ_DISAGREEMENT':concrete.length>1?'MULTI_READ_AGREEMENT':calls.length?'SINGLE_READ_SUPPORT':'NO_READ_SUPPORT'};
   });
   activeDisagreements=activeEvidence.filter(x=>String(x.Disagreement).toLowerCase()==='yes');
   // Navigation must include both concrete A/C/G/T disagreements and positions
   // where an ambiguous call (N or IUPAC) is paired with a concrete call.
   // The latter are review sites, not unambiguous nucleotide disagreements.
   activeReviewPositions=activeEvidence.filter(x=>{
     const observed=String(x.ObservedBases||'').split(';').filter(b=>b && b!=='-');
     return new Set(observed).size>1 && observed.length>1;
   });
   const seqNode=$('#consensusSeq');
   if(seqNode)seqNode.innerHTML=[...activeConsensus].map((b,i)=>`<span class="base ${b==='N'?'ambiguous':''} ${String(activeEvidence[i]?.Disagreement).toLowerCase()==='yes'?'discordant':''}" data-active-pos="${i}" title="Position ${i+1}; ${esc(activeEvidence[i]?.ObservedBases||'')}; ${esc(activeEvidence[i]?.EvidenceClass||'')}">${b}</span>`).join('');
   $('#readAlignmentPanel .panel-heading h3').textContent=`Multiple-read alignment · ${candidate.CandidateID}`;
   $('#readAlignmentPanel .panel-heading p').textContent=`Reads included in ${candidate.CandidateID} are projected onto its exact stored sequence. Disagreements are recalculated for this candidate.`;
   $('#readAlignmentPanel .panel-heading .small').textContent=`${alignmentReads.length} mapped reads`;
   $('#alignmentStart').max=Math.max(1,activeConsensus.length);
   $('#alignmentRange').max=Math.max(1,activeConsensus.length);
   $('#alignmentStatus').textContent=alignmentReads.length?`${alignmentReads.length} candidate read(s) projected. ${activeDisagreements.length} unambiguous disagreements; ${activeReviewPositions.length} total read-conflict/ambiguous review positions.`:'No candidate reads could be aligned to this sequence.';
   showAlignmentPosition(0);
   renderActiveSequence();
 }
 document.querySelectorAll('.candidate-select').forEach(btn=>btn.addEventListener('click',()=>{
   const candidate=candidateById.get(btn.dataset.candidate); if(!candidate||!candidate.Sequence)return;
   selectedCandidate=candidate;
   activateCandidateAlignment(candidate);
   document.querySelectorAll('.candidate-card').forEach(card=>card.classList.remove('chosen'));
   btn.closest('.candidate-card').classList.add('chosen');
   selectedStatus.textContent=`Selected ${candidate.CandidateID} · ${candidate.len.toLocaleString()} bp · ${candidate.QCClass}`;
   selectedDownload.disabled=false;
 }));
 selectedDownload.onclick=()=>{if(!selectedCandidate)return;const name=`${r.sequence_id}|${selectedCandidate.CandidateID}|${selectedCandidate.QCClass}`;dl(`${r.sequence_id}_${selectedCandidate.CandidateID}.fasta`,`>${name}\n${(selectedCandidate.Sequence.match(/.{1,70}/g)||[]).join('\n')}\n`,'text/x-fasta');selectedAssemblySequences[r.sequence_id]={candidate:selectedCandidate,sequence:selectedCandidate.Sequence};$('#downloadSelectedBatch').disabled=false;$('#downloadSelectedBatch').textContent=`Download selected isolates FASTA (${Object.keys(selectedAssemblySequences).length})`;};
 $('#downloadSelectedBatch').onclick=()=>{const records=Object.entries(selectedAssemblySequences).map(([id,v])=>`>${id}|${v.candidate.CandidateID}|${v.candidate.QCClass}\n${(v.sequence.match(/.{1,70}/g)||[]).join('\n')}`).join('\n');if(records)dl('SangerStudio_selected_assemblies.fasta',records+'\n','text/x-fasta');};
 renderActiveSequence();

 let alignmentReads=r.read_alignment||[];
 let activeConsensus=consensus;
 const alignmentWidth=80;
 let alignmentPosition=0;
 function renderActiveSequence(){
   const node=$('#consensusSeq'); if(!node)return;
   node.querySelectorAll('.base').forEach((el,i)=>el.addEventListener('click',()=>{
     const x=activeEvidence[i]||{};
     $('#baseInspector').innerHTML=`<b>Position ${i+1} · ${esc(x.ConsensusBase||'')}</b> &nbsp; Supporting reads: ${esc(x.SupportingReads||'—')} &nbsp; Observed: ${esc(x.ObservedBases||'—')} &nbsp; <span>${esc(x.EvidenceClass||'')}</span>`;
     showAlignmentPosition(i);
   }));
 }
 function showAlignmentPosition(position){
   if(!activeConsensus.length)return;
   alignmentPosition=Math.max(0,Math.min(activeConsensus.length-1,position));
   const start=Math.floor(alignmentPosition/alignmentWidth)*alignmentWidth;
   const end=Math.min(activeConsensus.length,start+alignmentWidth);
   const grid=$('#readAlignment'); if(!grid)return;
   const rows=[{read_id:'CONSENSUS',status:r.qc_class,aligned_sequence:activeConsensus,consensus:true},...alignmentReads];
   grid.innerHTML=`<div class="alignment-grid" style="--cols:${end-start}">${rows.map(row=>{
     const label=row.consensus?'Consensus':`${row.read_id} · ${row.primer||''} · ${row.direction||''} · ${row.status||''}`;
     const cells=[...String(row.aligned_sequence||'').slice(start,end)].map((base,j)=>{
       const pos=start+j, ref=activeConsensus[pos]||'N';
       const cls=base==='-'?'gap':base==='N'||!('ACGT'.includes(base))?'ambig':(!row.consensus&&'ACGT'.includes(ref)&&base!==ref?'diff':'');
       return `<span class="aln-cell ${cls} ${row.consensus?'aln-consensus':''} ${pos===alignmentPosition?'selected':''}" data-pos="${pos}" title="Position ${pos+1}: ${esc(row.read_id)} = ${esc(base)}">${esc(base)}</span>`;
     }).join('');
     return `<div class="aln-label">${esc(label)}</div>${cells}`;
   }).join('')}</div>`;
   $('#alignmentStart').value=alignmentPosition+1;
   $('#alignmentRange').value=alignmentPosition+1;
   $('#alignmentWindowLabel').textContent=`${start+1}–${end} / ${activeConsensus.length}`;
   grid.querySelectorAll('.aln-cell').forEach(cell=>cell.addEventListener('click',()=>{
     const pos=Number(cell.dataset.pos); alignmentPosition=pos;
     const x=activeEvidence[pos]||{};
     $('#baseInspector').innerHTML=`<b>Position ${pos+1} · ${esc(x.ConsensusBase||'')}</b> &nbsp; Supporting reads: ${esc(x.SupportingReads||'—')} &nbsp; Observed: ${esc(x.ObservedBases||'—')} &nbsp; <span>${esc(x.EvidenceClass||'')}</span>`;
     showAlignmentPosition(pos);
   }));
 }
 if($('#alignmentRange')){
   $('#alignmentRange').oninput=e=>showAlignmentPosition(Number(e.target.value)-1);
   $('#alignmentStart').onchange=e=>showAlignmentPosition(Number(e.target.value)-1);
   $('#nextDisagreement').onclick=()=>{
     if(!activeReviewPositions.length){$('#alignmentStatus').textContent='No positions with differing read calls or ambiguous-vs-concrete evidence were found in this candidate alignment. Candidate overlap-conflict counts are a separate assembly/QC metric.';return;}
     const next=activeReviewPositions.find(x=>Number(x.ConsensusPosition)>alignmentPosition+1);
     const target=next?Number(next.ConsensusPosition)-1:Number(activeReviewPositions[0].ConsensusPosition)-1;
     showAlignmentPosition(target);
     const x=activeEvidence[target]||{};
     $('#baseInspector').innerHTML=`<b>Position ${target+1} · ${esc(x.ConsensusBase||'')}</b> &nbsp; Supporting reads: ${esc(x.SupportingReads||'—')} &nbsp; Observed: ${esc(x.ObservedBases||'—')} &nbsp; <span>${esc(x.EvidenceClass||'')}</span>`;
   };
   showAlignmentPosition(0);
 }
 renderTaxonomyResults(r);
 $('#images').innerHTML=r.artifacts.map(a=>`<figure><a href="${a.url}" target="_blank"><img loading="lazy" src="${a.url}" alt="${esc(a.label)}"></a><figcaption>${esc(a.label)}</figcaption></figure>`).join('');$('#reportLinks').innerHTML=r.files.map(f=>`<a href="${f.url}" target="_blank">${esc(f.label)}</a>`).join('');
 let review=r.qc_class==='PASS_WITH_WARNINGS'&&r.recommendation==='MANUAL_REVIEW';$('#reviewBox').classList.toggle('hidden',!review);$('#reviewStatus').textContent='';if(!review&&r.qc_class==='FAIL_QC')$('#why').innerHTML+=`<br><br><b>This marker failed automated QC.</b> Its eligibility for any exception decision is controlled by the attached cohort review template; inspect the sequence evidence and use the cohort review panel where appropriate.`;$('#results').scrollIntoView({behavior:'smooth',block:'start'})
}
async function saveReview(decision){if(!lastInspection)return;try{let r=await api('/api/review',{sequence_id:lastInspection.sequence_id,decision,reviewer:$('#reviewer').value,notes:$('#reviewNotes').value});$('#reviewStatus').textContent=`Saved ${decision} decision to ${r.path}`;$('#reviewFile').href='/artifact?path=marker_review_decisions.tsv'}catch(e){$('#reviewStatus').textContent=e.message;$('#reviewStatus').className='error'}}
async function init(){try{studioState=await api('/api/state');$('#outputBadge').textContent='Output: '+studioState.out;$('#outputPath').textContent=studioState.out;$('#mapPath').value=studioState.map_path||'';$('#screenRefPath').value=studioState.screen_ref||'';$('#screenTaxaPath').value=studioState.screen_taxa||'';updateTaxonomyScreenStatus(studioState.taxonomy_screen_enabled?'Reference available for optional VSEARCH re-screening.':'No reference configured for optional re-screening.');renderOnboardingTaxonomyState(studioState.onboarding_taxonomy);renderSettings();renderSamples();$('#sequenceType').onchange=applySequenceProfile;applySequenceProfile();$('#chooseInputFolderBtn').onclick=chooseInputFolder;$('#chooseMapBtn').onclick=chooseMap;$('#chooseReadsBtn').onclick=chooseReads;$('#chooseOutputBtn').onclick=chooseOutput;$('#shutdownBtn').onclick=shutdownStudio;$('#chooseTaxonomyRefBtn').onclick=chooseTaxonomyRef;$('#chooseTaxonomyTableBtn').onclick=chooseTaxonomyTable;$('#applyTaxonomySettingsBtn').onclick=applyTaxonomySettings;$('#rerunTaxonomyCheck').onchange=()=>updateTaxonomyScreenStatus($('#rerunTaxonomyCheck').checked?'Fresh VSEARCH will run on the current trimmed reads during the next inspection.':'Off: use pre-screened onboarding taxonomy where available.');$('#attachReviewBtn').onclick=attachReview;$('#loadMap').onclick=()=>loadMap($('#mapPath').value);$('#sampleSearch').oninput=renderSamples;$('#runBtn').onclick=inspect;$('#runBtnSidebar').onclick=inspect;$('#processAllBtn').onclick=processAll;$('#resetSettings').onclick=renderSettings;$('#approveBtn').onclick=()=>saveReview('approve');$('#rejectBtn').onclick=()=>saveReview('reject');$('#cohortReviewSearch').oninput=renderCohortReviewRows;$('#cohortReviewFilter').onchange=renderCohortReviewRows;$('#cohortReviewRows').addEventListener('change',e=>{let tr=e.target.closest('tr[data-seq]');if(!tr||!cohortReview)return;let record=cohortReview.rows.find(r=>r.sequence_id===tr.dataset.seq);if(!record)return;if(e.target.classList.contains('cohortDecision')){record.decision=e.target.value;let n=cohortReview.rows.filter(r=>String(r.decision||'').trim()).length;$('#cohortReviewPill').textContent=`${n}/${cohortReview.total} decided`}if(e.target.classList.contains('cohortNotes'))record.notes=e.target.value;if(e.target.classList.contains('fallbackReadChoice')){record.selected_read_id=e.target.value;record.read_selection_source=e.target.value&&e.target.value!==record.recommended_read_id?'MANUAL':'AUTO_LONGEST_USABLE_READ';}});$('#cohortReviewRows').addEventListener('input',e=>{if(!e.target.classList.contains('cohortNotes')||!cohortReview)return;let tr=e.target.closest('tr[data-seq]'),record=cohortReview.rows.find(r=>r.sequence_id===tr?.dataset.seq);if(record)record.notes=e.target.value});$('#cohortReviewRows').addEventListener('click',e=>{let b=e.target.closest('[data-review-inspect]');if(!b)return;let id=b.dataset.reviewInspect;if((studioState.samples||[]).some(s=>s.sequence_id===id)){selectSample(id);$('#sampleHeading').scrollIntoView({behavior:'smooth',block:'start'})}else{setStatus($('#cohortReviewSaveStatus'),`${id} is not present in the currently loaded raw-read map. Attach the relevant partner sample map to inspect its chromatograms.`,true)}});$('#saveCohortReviewBtn').onclick=saveCohortReview;if(studioState.marker_review)renderCohortReview(studioState.marker_review);if(studioState.samples?.length){$('#mapStatus').textContent=studioState.onboarding_taxonomy?.loaded?studioState.onboarding_taxonomy.message:`Loaded ${studioState.samples.length} isolate(s) from the startup state.`;selectSample(studioState.samples[0].sequence_id)}}catch(e){$('#outputBadge').textContent='Studio connection failed';$('#mapStatus').textContent=e.message}}
init();
