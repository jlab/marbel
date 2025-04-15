from marbel.io_utils import is_bedtools_available, is_cat_available
from marbel.io_utils import concat_bed_files


def test_check_bedtools_present(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda _: "/usr/bin/bedtools")
    assert is_bedtools_available() is True


def test_check_bedtools_missing(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda _: None)
    assert is_bedtools_available() is False


def test_check_cat_present(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda _: "/usr/bin/cat")
    assert is_cat_available() is True


def test_check_cat_missing(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda _: None)
    assert is_cat_available() is False


def test_concat_bed_files(tmp_path):
    file1 = tmp_path / "file1.bed"
    file2 = tmp_path / "file2.bed"
    file1.write_text("chr1\t100\t200\n")
    file2.write_text("chr2\t300\t400\n")

    output_file = tmp_path / "output.bed"
    concat_bed_files(tmp_path, output_file)

    expected = "chr1\t100\t200\nchr2\t300\t400\n"
    assert output_file.read_text() == expected
