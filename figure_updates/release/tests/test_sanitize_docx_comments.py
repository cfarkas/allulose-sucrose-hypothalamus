#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path


RELEASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RELEASE))

import sanitize_docx_comments as sanitizer  # noqa: E402


CONTENT_TYPES = b"""<?xml version="1.0" encoding="UTF-8"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Override PartName="/word/document.xml" ContentType="document"/>
<Override PartName="/word/comments.xml" ContentType="comments"/>
<Override PartName="/word/people.xml" ContentType="people"/>
</Types>"""
RELATIONSHIPS = b"""<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://example/officeDocument"/>
<Relationship Id="rId2" Type="http://example/comments" Target="comments.xml"/>
<Relationship Id="rId3" Type="http://example/people" Target="people.xml"/>
</Relationships>"""
DOCUMENT = b"""<?xml version="1.0" encoding="UTF-8"?>
<w:document xmlns:w="urn:w"><w:body><w:p>
<w:commentRangeStart w:id="1"/><w:r><w:t>Scientific text</w:t></w:r>
<w:ins w:id="2" w:author="Private Reviewer"><w:r><w:t> accepted</w:t></w:r></w:ins>
<w:del w:id="3" w:author="Private Reviewer"><w:r><w:delText> deleted</w:delText></w:r></w:del>
<w:commentRangeEnd w:id="1"/><w:r><w:rPr/><w:commentReference w:id="1"/></w:r>
</w:p></w:body></w:document>"""
COMMENTS = b"""<?xml version="1.0" encoding="UTF-8"?>
<w:comments xmlns:w="urn:w"><w:comment w:id="1" w:author="Private Reviewer">
<w:p><w:r><w:t>private@example.org</w:t></w:r></w:p></w:comment></w:comments>"""
PEOPLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<w15:people xmlns:w15="urn:w15"><w15:person w15:author="Private Reviewer"
w15:userId="private-account-uuid"/></w15:people>"""
STYLES = b"""<?xml version="1.0" encoding="UTF-8"?>
<w:styles xmlns:w="urn:w"><w:style w:styleId="Normal"/>
<w:style w:styleId="CommentText1"><w:name w:val="Comment Text1"/></w:style>
</w:styles>"""


class OfficeSanitizerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="office-sanitize-test-")
        self.base = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def write_zip(path: Path, members: dict[str, bytes]) -> None:
        with zipfile.ZipFile(path, "x", compression=zipfile.ZIP_DEFLATED) as archive:
            for name, payload in members.items():
                archive.writestr(name, payload)

    def test_docx_comment_parts_markers_relationships_and_identities_are_removed(self):
        source = self.base / "reviewed.docx"
        output = self.base / "clean.docx"
        self.write_zip(
            source,
            {
                "[Content_Types].xml": CONTENT_TYPES,
                "word/document.xml": DOCUMENT,
                "word/_rels/document.xml.rels": RELATIONSHIPS,
                "word/comments.xml": COMMENTS,
                "word/people.xml": PEOPLE,
                "word/styles.xml": STYLES,
                "word/media/image.bin": b"unchanged",
            },
        )
        report = sanitizer.sanitize_package(source, output)
        with zipfile.ZipFile(output) as archive:
            names = archive.namelist()
            joined = b"\n".join(archive.read(name) for name in names)
            self.assertNotIn("word/comments.xml", names)
            self.assertNotIn("word/people.xml", names)
            self.assertNotIn(b"commentRange", joined)
            self.assertNotIn(b"commentReference", joined)
            self.assertNotIn(b"Private Reviewer", joined)
            self.assertNotIn(b"private@example.org", joined)
            self.assertIn(b"Scientific text", joined)
            self.assertIn(b" accepted", joined)
            self.assertNotIn(b" deleted", joined)
            self.assertEqual(archive.read("word/media/image.bin"), b"unchanged")
        self.assertEqual(report["kind"], "docx")
        self.assertEqual(len(report["removed_parts"]), 2)

    def test_xlsx_private_metadata_is_removed_but_sheet_content_is_exact(self):
        source = self.base / "private.xlsx"
        output = self.base / "clean.xlsx"
        workbook = b"""<?xml version="1.0" encoding="UTF-8"?>
<workbook xmlns="urn:s" xmlns:ac="urn:abs" xmlns:xr="urn:revision">
<ac:absPath url="https://onedrive.example/account/private/path"/>
<xr:revisionPtr documentId="doc-uuid" xr:uidLastSave="user-uuid"/>
<bookViews><workbookView xr:uid="view-uuid"/></bookViews>
</workbook>"""
        sheet = b"""<?xml version="1.0" encoding="UTF-8"?>
<worksheet xmlns="urn:s" xmlns:xr="urn:revision" xr:uid="sheet-uuid">
<sheetData><row r="1"><c r="A1" t="inlineStr"><is><t>42</t></is></c></row></sheetData>
</worksheet>"""
        core = b"""<?xml version="1.0" encoding="UTF-8"?>
<cp:coreProperties xmlns:cp="urn:cp"><cp:creator>Author</cp:creator>
<cp:lastModifiedBy>Private Editor</cp:lastModifiedBy></cp:coreProperties>"""
        self.write_zip(
            source,
            {
                "[Content_Types].xml": b"<Types/>",
                "xl/workbook.xml": workbook,
                "xl/worksheets/sheet1.xml": sheet,
                "docProps/core.xml": core,
                "xl/styles.xml": b"<styleSheet/>",
            },
        )
        sanitizer.sanitize_package(source, output)
        with zipfile.ZipFile(output) as archive:
            joined = b"\n".join(archive.read(name) for name in archive.namelist())
            self.assertNotIn(b"absPath", joined)
            self.assertNotIn(b"revisionPtr", joined)
            self.assertNotIn(b"lastModifiedBy", joined)
            self.assertNotIn(b"sheet-uuid", joined)
            self.assertNotIn(b"Private Editor", joined)
            self.assertIn(b">42<", joined)

    def test_in_place_batch_creates_byte_identical_hierarchical_backup(self):
        root = self.base / "Paper"
        source = root / "nested" / "reviewed.docx"
        source.parent.mkdir(parents=True)
        self.write_zip(
            source,
            {
                "[Content_Types].xml": CONTENT_TYPES,
                "word/document.xml": DOCUMENT,
                "word/_rels/document.xml.rels": RELATIONSHIPS,
                "word/comments.xml": COMMENTS,
            },
        )
        original = source.read_bytes()
        reports = sanitizer.sanitize_in_place(
            root,
            [Path("nested/reviewed.docx")],
            Path("recovery_support/backups"),
        )
        backup = root / "recovery_support/backups/nested/reviewed.docx"
        self.assertEqual(backup.read_bytes(), original)
        self.assertNotEqual(source.read_bytes(), original)
        manifest = json.loads(
            (root / "recovery_support/backups/backup_manifest.json").read_text()
        )
        self.assertEqual(manifest["files"][0]["path"], "nested/reviewed.docx")
        self.assertEqual(len(reports), 1)


if __name__ == "__main__":
    unittest.main()
