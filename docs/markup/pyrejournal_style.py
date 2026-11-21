# -*- python -*-
# -*- coding: utf-8 -*-
#
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2026 The PyLith authors and contributors
# SPDX-License-Identifier: MIT


"""Pygments style for PyLith output, matching the terminal color scheme."""

from pygments.style import Style
from pygments.token import Token


class PyreJournalStyle(Style):
    background_color = "#aaaaaa"
    styles = {
        Token.Generic.Marker: "#aaaaaa",  # >>
        Token.Keyword: "bold #00aa00",  # info (green)
        Token.Generic.Warning: "bold #ff6600",  # warning (orange)
        Token.Generic.Error: "bold #aa0000",  # error (red)
        Token.Name.Tag: "bold #aa00aa",  # (category) (purple)
        Token.Comment: "#aaaaaa",  # -- / dim text
        Token.Comment.Single: "italic #6666ff",  # # comments (inline/full)
        Token.Generic.Path: "#bbbbbb",  # paths / signatures
        Token.Text: "#000000",  # default text
    }


# End of file
