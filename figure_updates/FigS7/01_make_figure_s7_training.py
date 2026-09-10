#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
for parent in HERE.parents:
 candidate=parent/'scripts/shared'
 if (candidate/'render_method_figures.py').exists():
  sys.path.insert(0,str(candidate));paper_root=parent;break
from render_method_figures import main
if __name__=='__main__':main('S7',paper_root/'FigS7')
