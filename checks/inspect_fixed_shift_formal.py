"""Inspect the public Prove2Me pages. Optional local HTML cache avoids network access."""
import argparse
from hashlib import sha256
from html.parser import HTMLParser
from pathlib import Path
from urllib.request import Request, urlopen


NODES = (
    ("top", "de5a0a3d-3615-5c54-b970-82a08423d350"),
    ("node-premise", "ecb0f475-7d39-4e1d-bf96-947e13e6db6a"),
    ("node-p41", "298ec1ce-a5c8-4485-9da8-d1b541149a62"),
    ("node-l42", "d10b17f2-9d2a-458a-bb30-e9bd5b0834aa"),
    ("node-t43", "e42e2812-6311-4551-8c39-6a3d880cf5be"),
    ("node-t43b", "6c9f5dd7-cb2d-4b0d-9451-1ce69bb2040a"),
    ("node-l52", "5c309fe4-3a1c-48c0-972a-bf272ca10df9"),
)


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.field = None
        self.end = None
        self.blocks = []
        self.info = {}

    def handle_starttag(self, tag, attrs):
        css = dict(attrs).get("class", "").split()
        field = next((k for k in ("h-pedia-leanname", "h-pedia-status") if k in css), None)
        if field or tag == "pre":
            self.field = field or "pre"
            self.end = tag
            if self.field == "pre":
                self.blocks.append("")

    def handle_endtag(self, tag):
        if tag == self.end:
            self.field = self.end = None

    def handle_data(self, data):
        if self.field == "pre":
            self.blocks[-1] += data
        elif self.field:
            self.info[self.field] = self.info.get(self.field, "") + data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, help="directory containing lzr-<key>.html")
    args = parser.parse_args()
    for key, uid in NODES:
        if args.cache:
            source = (args.cache / f"lzr-{key}.html").read_text(encoding="utf-8-sig")
        else:
            req = Request(f"https://prove2.me/theorems/{uid}", headers={"User-Agent": "Research-audit/1"})
            with urlopen(req, timeout=30) as response:
                source = response.read().decode("utf-8")
        page = Page()
        page.feed(source)
        assert page.blocks, f"Missing public Lean block at {uid}"
        lean = "\n".join(block.replace("\r\n", "\n").strip() for block in page.blocks)
        if key == "node-premise":
            for field in ("closed_schedule_contraction", "one_step_packet_hermite_mode"):
                assert field in lean
        print(uid, page.info.get("h-pedia-status", "?"), page.info.get("h-pedia-leanname", "?"),
              "public-code-sha256", sha256(lean.encode()).hexdigest())
    print("Hashes cover displayed statements/definitions. No accepted proof replay or transitive axiom inventory.")


if __name__ == "__main__":
    main()
