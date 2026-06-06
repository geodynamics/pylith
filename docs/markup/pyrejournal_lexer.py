"""Pygments lexer for PyLith simulation output (semantic, prefix-driven)."""

from pygments.lexer import RegexLexer, bygroups
from pygments.token import Token

# Semantic tokens
Marker  = Token.Generic.Marker       # >>
Info    = Token.Keyword              # info (green)
Warn    = Token.Generic.Warning      # warning (orange)
Err     = Token.Generic.Error        # error (red)
Channel = Token.Name.Tag             # (application-flow) etc. (magenta)
Dash    = Token.Comment              # -- prefix and dim text
Path    = Token.Generic.Path         # file paths / C++ signatures
Comment = Token.Comment.Single       # # ... comments (inline or full-line)
Proc    = Token.Generic.Path         # [Proc]


class PyreJournalLexer(RegexLexer):
    name = 'pyrejournal'
    aliases = ['pyrejournal', "pylith_output"]

    tokens = {
        # Reusable: a '#' (and the rest of the line) is a comment. Placed at the
        # end of body-bearing states so prefix rules win first.
        'comment': [
            (r'#.*?$', Comment),
        ],
        'root': [
            # >> source/file reference lines (paths don't take # comments)
            (r'^\s*(>>)( )(.*)$', bygroups(Marker, Token.Text, Path)),
            # -- <level> (<category>) — level colored, category magenta.
            # Stop the match after the category so any trailing text (incl. an
            # inline # comment) is handled by the 'body' state.
            (r'^\s*(-- )(info)( )(\([\w._\-]+\))',
                bygroups(Dash, Info, Token.Text, Channel), 'body'),
            (r'^\s*(-- )(warning)( )(\([\w._\-]+\))',
                bygroups(Dash, Warn, Token.Text, Channel), 'body'),
            (r'^\s*(-- )(error)( )(\([\w._\-]+\))',
                bygroups(Dash, Err, Token.Text, Channel), 'body'),
            # plain -- lines: emit the prefix, then let 'body' handle the rest
            (r'^\s*(--)( )', bygroups(Dash, Token.Text), 'body'),
            (r'^\s*(--)$', Dash),
            # PETSC ERROR
            #(r'^\s*(\[[0-9]+\])(PETSC ERROR: #[0-9]+.*)$',
            #    bygroups(Proc, Path)),
            #(r'^\s*(\[[0-9]+\])(PETSC ERROR:)( )',
            #    bygroups(Proc, Err, Token.Text), 'body'),
            (r'^\s*(\[[0-9]+\])(PETSC ERROR: #[0-9]+.*)$',
                bygroups(Token.Text, Token.Text)),
            (r'^\s*(\[[0-9]+\])(PETSC ERROR:)( )',
                bygroups(Token.Text, Token.Text, Token.Text), 'body'),
            # any other line: hand off to 'body' so inline # comments are caught
            (r'', Token.Text, 'body'),
        ],
        # Body of a line: text until a '#', then the comment, until newline.
        'body': [
            (r'\n', Token.Text, '#pop'),
            (r'#.*?(?=\n|$)', Comment),
            (r'[^#\n]+', Token.Text),
        ],
    }
