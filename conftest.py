"""
The test suite runs from the source tree against an *editable* install
(`pip install -e .`), as upstream orvara always has: the test configs address
their data by paths relative to the repository root, and `orvara/tests` is a
package inside the source tree.  With a non-editable install the source tree
shadows the installed package and has no compiled `orvara.orbit`, which
produces a baffling ImportError -- turn it into a clear message instead.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))


def pytest_sessionstart(session):
    src = os.path.join(_HERE, 'orvara')
    if not any(name.startswith('orbit.') and name.endswith(('.so', '.pyd'))
               for name in os.listdir(src)):
        sys.stderr.write(
            '\n*** orvara/orbit.* is not compiled in this source tree. The test '
            'suite runs against an editable install:\n'
            '***     pip install -e ".[test]"\n'
            '*** (a plain `pip install .` puts the extension in site-packages, '
            'where these tests cannot see it).\n\n')
