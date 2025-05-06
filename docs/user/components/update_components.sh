~/src/cig/pythia/bin/pyre_doc_components.py --package=PyLith --src-path=$PYENV/pylith-debug/lib/python3.12/site-packages/pylith --output-path=./ --manual-toc=toc.json
cat index.md | sed -e "s#apps/index.md#implementations.md\\napps/index.md#g" > tmp.md && mv -f tmp.md index.md

# Must edit PyLithApp.md and ConvertMeshApp.md manually.