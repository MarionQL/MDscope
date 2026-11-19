"""Launch the optional authoring app on localhost."""
import sys
from pathlib import Path


def main():
    try:
        from streamlit.web import cli
    except ImportError:
        raise SystemExit('Install Studio with: pip install "mdscope[studio]"') from None
    sys.argv = ["streamlit", "run", str(Path(__file__).with_name("app.py")),
                "--server.address=127.0.0.1", *sys.argv[1:]]
    return cli.main()
