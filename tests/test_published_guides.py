"""Keep release identifiers validated after shortening the landing page."""
import hashlib
from pathlib import Path
import shutil
import tempfile
import unittest

from repository_tools import PublicationError, validate_repository


class PublishedGuideTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        source = Path(__file__).resolve().parents[1]
        for pattern in ("*.md", "*.py", "*.txt", "*.json", "*.env", "LICENSE", "LICENSES/*", ".github/workflows/*.yml"):
            for path in source.glob(pattern):
                if not path.is_file():
                    continue
                target = self.root / path.relative_to(source)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, target)

    def test_beginner_page_and_linked_details_preserve_release_contract(self):
        self.assertEqual(validate_repository(self.root)["archive_count"], 10)

    def test_missing_detailed_checksum_still_fails_publication(self):
        guide = self.root / "DETAILED_GUIDE.md"
        digest = hashlib.sha256((self.root / "manifest.json").read_bytes()).hexdigest()
        guide.write_text(guide.read_text().replace(digest, "omitted digest"))
        with self.assertRaisesRegex(PublicationError, "omit required resolved values"):
            validate_repository(self.root)

    def test_unlinked_details_cannot_hide_required_release_information(self):
        readme = self.root / "README.md"
        readme.write_text(readme.read_text().replace("(DETAILED_GUIDE.md)", "(unavailable.md)"))
        with self.assertRaisesRegex(PublicationError, "must link"):
            validate_repository(self.root)

    def test_missing_beginner_launcher_blocks_publication(self):
        (self.root / "reproduce.py").unlink()
        with self.assertRaisesRegex(PublicationError, "Required beginner-guide file is absent"):
            validate_repository(self.root)


if __name__ == "__main__":
    unittest.main()
