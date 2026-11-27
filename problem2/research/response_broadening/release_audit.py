"""Read-only GitHub Release evidence and historical NPZ comparison."""
import hashlib,io,json,subprocess,zipfile
from pathlib import Path
import numpy as np
from ...src.delay_audit import write_json,write_csv,compare_tree

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=Path.home()/'.cache/apmcm-b/releases'


def inspect(out,refresh=False):
    inventory=out/'release_inventory.json'
    if refresh or not inventory.exists():
        releases=[]
        for tag in ['problem1-final-data-v1','problem2-final-data']:
            releases.append(json.loads(subprocess.check_output(['gh','release','view',tag,'--repo','ma39nong/Steelmaking-process','--json','tagName,name,body,url,assets'],text=True)))
        write_json(inventory,releases)
    else: releases=json.loads(inventory.read_text())
    package=CACHE/'problem2-final-data/Steelmaking-process-q2-final-data.zip'
    checksum=CACHE/'problem2-final-data/Steelmaking-process-q2-final-data.SHA256'
    q1=CACHE/'problem1-final-data-v1/data_package_manifest.json'
    sums=CACHE/'problem1-final-data-v1/SHA256SUMS'
    if not all(p.exists() for p in [package,checksum,q1,sums]):
        raise FileNotFoundError('Download named Release checksum/manifest assets to '+str(CACHE)+' first; see README')
    digest=lambda b:hashlib.sha256(b).hexdigest()
    expected=checksum.read_text().split()[0];actual=digest(package.read_bytes())
    if expected!=actual:raise AssertionError('Release package checksum mismatch')
    q1lines=[line.split() for line in sums.read_text().splitlines() if line.strip()]
    q1sha=next(words[0] for words in q1lines if words[-1].lstrip('*')=='data_package_manifest.json')
    if digest(q1.read_bytes())!=q1sha:raise AssertionError('first-question manifest checksum mismatch')
    q1manifest=json.loads(q1.read_text());rows=[]
    with zipfile.ZipFile(package) as z:
        names=z.namelist();checks=z.read('SHA256SUMS').decode()
        for line in checks.splitlines():
            fields=line.split(maxsplit=1)
            if len(fields)!=2:continue
            filename=fields[1].strip().lstrip('*').removeprefix('./')
            if digest(z.read(filename))!=fields[0]:raise AssertionError('member checksum mismatch: '+filename)
        for name in names:
            if not name.endswith('.npz'):continue
            with np.load(io.BytesIO(z.read(name)),allow_pickle=False) as old,np.load(ROOT/'problem2/final/delay_audit/baseline/data'/Path(name).name,allow_pickle=False) as new:
                if set(old.files)!=set(new.files):raise AssertionError('historical NPZ field set mismatch')
                for field in old.files:
                    a,b=old[field],new[field];same=a.shape==b.shape
                    numeric=a.dtype.kind in 'biufc'
                    exact=bool(np.array_equal(a,b,equal_nan=True)) if numeric else bool(np.array_equal(a,b))
                    near=bool(same and np.allclose(a,b,atol=1e-8,rtol=1e-8,equal_nan=True)) if numeric else exact
                    if field=='metadata':near=not compare_tree(json.loads(str(a)),json.loads(str(b)))
                    maximum=float(np.max(abs(a[np.isfinite(a)]-b[np.isfinite(a)]))) if numeric and same and np.isfinite(a).any() and np.array_equal(np.isfinite(a),np.isfinite(b)) and a.dtype.kind!='b' else None
                    rows.append(dict(heat_id=Path(name).name.split('_')[0],field=field,shape_match=same,exact_match=exact,
                                     numerical_or_metadata_match=near,maximum_finite_difference=maximum))
    write_csv(out/'release_npz_comparison.csv',rows)
    protected_identity_fields=lambda r:('mask' in r['field'] or r['field'] in ['sequence','time'] or r['field'].endswith(('center_index','window_time')))
    if any(not r['numerical_or_metadata_match'] for r in rows) or any(not r['exact_match'] for r in rows if protected_identity_fields(r)):
        raise AssertionError('Release numerical/interface differences need investigation')
    write_json(out/'release_identity.json',dict(releases=[dict(tag=r['tagName'],url=r['url'],assets=[a['name'] for a in r['assets']]) for r in releases],
        q2_zip_sha256=actual,zip_checksum_verified=True,all_zip_member_checksums_verified=True,archive_inventory=names,
        historical_npz_heats=len({r['heat_id'] for r in rows}),historical_npz_fields=len(rows),
        historical_numeric_and_metadata_match=all(r['numerical_or_metadata_match'] for r in rows),
        masks_centers_and_time_axes_exact=True,
        q1_manifest_sha256=q1sha,q1_manifest_checksum_verified=True,
        q1_manifest_declares_raw_xlsb_included=q1manifest['original_xlsb_and_process_columns_included'],
        raw_xlsb_assets_listed=any('.xlsb' in a['name'].lower() for r in releases for a in r['assets']),
        designated_new_raw_attachment1_verified=False,
        distinction='Historical processed NPZ identity is now checked directly. Releases do not supply the designated new raw attachment1; these are separate checks.'))
    print('Release NPZ comparison verified',len(rows),'fields',flush=True)


if __name__=='__main__':
    out=HERE/'results';out.mkdir(exist_ok=True);inspect(out,refresh=True)
