"""Package frozen arrays; verify provenance without fitting or changing any model."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[3]
FINAL = ROOT / 'problem1/final'
OUTPUT = FINAL / 'delivery'
DATA_COMMIT = '3366db522887b2ca0f5f461ea6bd52ed71def1ba'


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def put_file(archive, path, name):
    # Fixed metadata: identical input files yield an identical archive.
    info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_STORED
    info.external_attr = 0o100644 << 16
    with archive.open(info, 'w', force_zip64=True) as dest, path.open('rb') as src:
        for block in iter(lambda: src.read(1024 * 1024), b''):
            dest.write(block)


def main():
    with (FINAL / 'results/artifact_manifest.csv').open() as stream:
        expected = {(r['heat_id'], r['artifact']): r for r in csv.DictReader(stream)}
    with (FINAL / 'results/final_summary.csv').open() as stream:
        summary = list(csv.DictReader(stream))
    if len(summary) != 10 or sum(int(r['repair_points']) for r in summary) != 1043:
        raise ValueError('Frozen ten-heat/1043-repair summary mismatch')
    payload = []
    for row in summary:
        heat = row['heat_id']
        for kind in ('interface', 'difference'):
            path = FINAL / 'data' / f'{heat}_{kind}.npz'
            expected_file = expected[(heat, kind)]
            digest = sha256(path)
            if digest != expected_file['sha256'] or path.stat().st_size != int(expected_file['bytes']):
                raise ValueError(f'Frozen artifact checksum mismatch: {path.name}')
            payload.append((path, digest))
    paths = [
        ROOT / 'data_manifest.csv', ROOT / 'problem1/v3/src/interface.py',
        FINAL / 'reports/数据交接与后续问题使用说明.md',
        OUTPUT / 'README.md', OUTPUT / 'read_example.py',
    ] + [FINAL / 'results' / name for name in (
        'final_summary.csv', 'actual_repairs.csv', 'statistical_repairs_evidence.csv',
        'native_mask_summary.csv', 'conditional_tiles.csv', 'artifact_manifest.csv',
    )]
    payload.extend((path, sha256(path)) for path in paths)
    entries = [{'path': path.relative_to(ROOT).as_posix(), 'bytes': path.stat().st_size,
                'sha256': digest} for path, digest in payload]
    manifest = {
        'dataset': 'Problem 1 frozen final recorded-scale offline spectra',
        'data_source_commit': DATA_COMMIT,
        'packaging_base_commit': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'heats': 10, 'actual_repairs': 1043,
        'coordinate_units': {'time': 'relative seconds', 'wavelength': 'nm', 'intensity': 'counts'},
        'offline': True, 'original_xlsb_and_process_columns_included': False,
        'archive': 'problem1-final-data.zip',
        'payload_bytes': sum(row['bytes'] for row in entries),
        'files': entries,
    }
    manifest_path = OUTPUT / 'data_package_manifest.json'
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    target = OUTPUT / manifest['archive']
    pending = OUTPUT / (target.name + '.partial')
    try:
        with zipfile.ZipFile(pending, 'w', allowZip64=True) as archive:
            for path, digest in payload:
                put_file(archive, path, path.relative_to(ROOT).as_posix())
            put_file(archive, manifest_path, manifest_path.relative_to(ROOT).as_posix())
        if pending.stat().st_size >= 2 * 1024**3:
            raise ValueError('Archive must be smaller than 2 GiB per Release asset')
        with zipfile.ZipFile(pending) as archive:
            for entry in entries:
                digest = hashlib.sha256()
                with archive.open(entry['path']) as stream:
                    for block in iter(lambda: stream.read(1024 * 1024), b''):
                        digest.update(block)
                if digest.hexdigest() != entry['sha256']:
                    raise ValueError(f'Archive checksum mismatch: {entry["path"]}')
            # Reading all entries above validates CRC as well as our SHA256.
            if archive.read(manifest_path.relative_to(ROOT).as_posix()) != manifest_path.read_bytes():
                raise ValueError('Embedded manifest mismatch')
        pending.replace(target)
    finally:
        pending.unlink(missing_ok=True)
    checksum_path = OUTPUT / 'SHA256SUMS'
    checksum_path.write_text(''.join(f'{sha256(p)}  {p.name}\n' for p in (target, manifest_path)))
    print(json.dumps({'archive': str(target), 'bytes': target.stat().st_size,
                      'MiB': round(target.stat().st_size / 1024**2, 2),
                      'payload_files': len(entries), 'arrays_verified': 20,
                      'archive_crc_and_sha256': 'passed', 'published_to_github': False},
                     ensure_ascii=False))


if __name__ == '__main__':
    main()
