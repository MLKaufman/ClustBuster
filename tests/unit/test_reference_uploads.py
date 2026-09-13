from pathlib import Path

import pandas as pd
import pytest

from clustbuster.resources.references.details import reference_compatibility
from clustbuster.services.imports import SessionFiles, UploadError
from clustbuster.services.reference_uploads import import_reference_upload


def test_uploaded_csv_with_rownames_becomes_session_reference(tmp_path: Path) -> None:
    source = tmp_path / 'input.csv'
    source.write_text(',T cell,B cell\nCD3D,4,0\nMS4A1,0,5\nNA,1,1\n')
    session = SessionFiles.create(tmp_path / 'sessions')
    loaded = import_reference_upload(
        {'name': '../../my-reference.csv', 'datapath': str(source)}, session,
        max_upload_mb=1, name='My atlas', species='Homo sapiens', normalization='log1p',
    )
    assert loaded.matrix_path.is_relative_to(session.uploads)
    assert loaded.summary.name == 'My atlas'
    assert loaded.summary.cell_types == ('T cell', 'B cell')
    assert loaded.summary.reference_id.startswith('upload:')
    assert loaded.metadata['normalization'] == 'log1p'
    assert len(loaded.checksum) == 64
    table = pd.read_csv(loaded.matrix_path, sep='\t', keep_default_na=False)
    assert table['gene'].tolist() == ['CD3D', 'MS4A1', 'NA']
    matches = reference_compatibility(loaded, ('CD3D', 'MS4A1', 'NA'))
    assert len(matches) == 3
    assert (matches['Match'] == 'Shared').all()
    assert source.read_text().startswith(',T cell,B cell')
    session.cleanup()
    assert not loaded.matrix_path.exists()
    assert source.exists()


@pytest.mark.parametrize('text, message', [
    ('gene\tT\tT\nCD3D\t1\t2\n', 'column names'),
    ('gene\tT\nCD3D\t1\nCD3D\t2\n', 'gene identifiers'),
    ('gene\tT\n\t1\n', 'gene identifiers'),
    ('gene\tT\nCD3D\t\n', 'numeric and finite'),
    ('gene\tT\nCD3D\tinf\n', 'numeric and finite'),
    ('gene\tT\nCD3D\t1\t2\n', 'different number of columns'),
    ('gene\tT\n', 'gene identifiers'),
])
def test_invalid_reference_is_rejected_and_cleaned(
    tmp_path: Path, text: str, message: str,
) -> None:
    source = tmp_path / 'input.tsv'
    source.write_text(text)
    session = SessionFiles.create(tmp_path / 'sessions')
    with pytest.raises(UploadError, match=message):
        import_reference_upload({'name': 'input.tsv', 'datapath': str(source)}, session,
                                max_upload_mb=1)
    assert not list(session.uploads.iterdir())


def test_upload_limit_and_format_are_enforced(tmp_path: Path) -> None:
    source = tmp_path / 'file'
    source.write_text('gene\tT\nCD3D\t1\n')
    session = SessionFiles.create(tmp_path / 'sessions')
    with pytest.raises(UploadError, match='CSV'):
        import_reference_upload({'name': 'matrix.rds', 'datapath': str(source)}, session,
                                max_upload_mb=1)
    with pytest.raises(UploadError, match='upload limit'):
        import_reference_upload({'name': 'matrix.tsv', 'datapath': str(source)}, session,
                                max_upload_mb=0)
