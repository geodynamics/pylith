% -*- LaTeX -*-
%
% ----------------------------------------------------------------------
%
% Brad T. Aagaard, U.S. Geological Survey
% Charles A. Williams, GNS Science
% Matthew G. Knepley, University of Chicago
%
% This code was developed as part of the Computational Infrastructure
% for Geodynamics (http://geodynamics.org).
%
% Copyright (c) 2010-2017 University of California, Davis
%
% See COPYING for license information.
%
% ----------------------------------------------------------------------
%
% pylithdoc.cls    Style file for PyLith documentation.
\NeedsTeXFormat{LaTeX2e}[1994/06/01]
\ProvidesClass{pylithdoc}[2017/02/17 v0.1 PyLith Documentation class]

\ProcessOptions*
\LoadClass[10pt,twoside]{book}

% Packages
\RequirePackage{fourier} % texlive-fonts-extra
\RequirePackage{float}
\RequirePackage{varioref} % \vref
\RequirePackage{color}
\RequirePackage{listings}
\RequirePackage{pifont}
\RequirePackage{calc}
\RequirePackage{tikz}

\RequirePackage[%
  unicode=true,%
  pdfusetitle,%
  bookmarks=true,
  bookmarksnumbered=false,%
  bookmarksopen=false,%
  breaklinks=true,%
  pdfborder={0 0 0},%
  backref=false,%
  colorlinks=true]{hyperref}

% ------------------------------------------------------------------
% Basic page layout
\setlength{\textheight}{9.0in}
\setlength{\textwidth}{7.0in}
\setlength{\topmargin}{-20pt}
\setlength{\headheight}{16pt}
\setlength{\headsep}{4pt}
\setlength{\oddsidemargin}{-0.25in}
\setlength{\evensidemargin}{-0.25in}
\setlength{\footskip}{24pt}
%
\setlength{\floatsep}{12pt}
\setlength{\textfloatsep}{12pt}
\setlength{\intextsep}{12pt}
\setlength{\abovecaptionskip}{0pt}
\setlength{\belowcaptionskip}{6pt}
%
\setlength{\parindent}{0.25in}
\setlength{\parskip}{6pt}
%
% ------------------------------------------------------------------
% Font Sizes
\RequirePackage{times}
\renewcommand{\LARGE}{\@setfontsize\LARGE{17}{22}}
\renewcommand{\Large}{\@setfontsize\Large{14.4}{18}}
\renewcommand{\large}{\@setfontsize\large{10.95}{13.6}}
\renewcommand{\normalsize}{\@setfontsize\normalsize{10}{12}}
\renewcommand{\small}{\@setfontsize\small{9}{10}}
\renewcommand{\footnotesize}{\@setfontsize\footnotesize{8}{9.5}}
%
% Colors
\definecolor{yellow}{rgb}{1.0, 1.0, 0.45} % 255/255/115
\definecolor{dkyellow}{rgb}{0.9, 0.9, 0.0} % % 230/230/0

\definecolor{ltorange}{rgb}{1.0, 0.74, 0.41} % 255/188/105
\definecolor{orange}{rgb}{0.96, 0.50, 0.0} % 246/127/0

\definecolor{ltred}{rgb}{1.0, 0.25, 0.25} % 255/64/64
\definecolor{red}{rgb}{0.79, 0.00, 0.01} % 201/0/3

\definecolor{ltpurple}{rgb}{0.81, 0.57, 1.00} % 206/145/255
\definecolor{purple}{rgb}{0.38, 0.00, 0.68} % 97/1/175

\definecolor{ltblue}{rgb}{0.2, 0.73, 1.0} % 51/187/255
\definecolor{blue}{rgb}{0.12, 0.43, 0.59} % 30/110/150

\definecolor{ltltgreen}{rgb}{0.7, 1.00, 0.7} % 96/204/14
\definecolor{ltgreen}{rgb}{0.37, 0.80, 0.05} % 96/204/14
\definecolor{green}{rgb}{0.23, 0.49, 0.03} % 59/125/8
  
\definecolor{dkslate}{rgb}{0.18, 0.21, 0.28} % 47/53/72
\definecolor{mdslate}{rgb}{0.45, 0.50, 0.68} % 114/127/173
\definecolor{ltslate}{rgb}{0.85, 0.88, 0.95} % 216/225/229

% ------------------------------------------------------------------
% Headings
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}

% ------------------------------------------------------------------
% Macros
\usetikzlibrary{shadows.blur}
% Warning
\newcommand{\warning}[1]{%
  \begin{center}\begin{tikzpicture}
      \node[draw=none,
        shade,
        text width=\textwidth-2.0in,
        inner sep=1em,
        top color=red!50,
        bottom color=red!5,
        rounded corners=1em,
        blur shadow={shadow blur steps=5}
      ]{\raisebox{1em}{{\LARGE\danger} {\bf Warning}}\\ #1};
    \end{tikzpicture}\end{center}}
%
% Tip
\newcommand{\tip}[1]{%
  \begin{center}\begin{tikzpicture}
      \node[draw=none,
        shade,
        text width=\textwidth-2.0in,
        inner sep=1em,
        top color=green!50,
        bottom color=green!5,
        rounded corners=1em,
        blur shadow={shadow blur steps=5}
      ]{\raisebox{1em}{{\LARGE\ding{72}} {\bf Tip}}\\ #1};
    \end{tikzpicture}\end{center}}
%
% Important
\newcommand{\important}[1]{%
  \begin{center}\begin{tikzpicture}
      \node[draw=none,
        shade,
        text width=\textwidth-2.0in,
        inner sep=1em,
        top color=blue!50,
        bottom color=blue!5,
        rounded corners=1em,
        blur shadow={shadow blur steps=5}
      ]{\raisebox{1em}{{\LARGE\ding{234}} {\bf Important}}\\ #1};
    \end{tikzpicture}\end{center}}

% New feature
\newcommand{\newfeature}[1]{%
  \begin{tikzpicture}
      \node[draw=none,
        inner sep=0.3em,
        top color=purple!60,
        bottom color=purple!20,
        rounded corners=0.2em,
        blur shadow={shadow blur steps=5}
      ]{\bf New in #1};
    \end{tikzpicture}}
%
% ------------------------------------------------------------------
% Code blocks
% .cfg
\lstdefinelanguage{cfg}{%
  basicstyle=\small\ttfamily,
  frame=tb,%
  morecomment=[l][\color{orange}\itshape]{\#},
  morecomment=[l][\color{orange}\itshape]{;},
  moredelim=[is][keywordstyle]{<h>}{</h>},
  moredelim=[is][\color{purple}\bfseries]{<p>}{</p>},
  moredelim=[is][\color{blue}\bfseries]{<f>}{</f>},
}
%
% shell  
\lstnewenvironment{shell}{\lstset{%
    language=bash,%
    frame=single,%
    basicstyle=\small\ttfamily,%
    commentstyle=\itshape\color{orange},
}}{}
\lstnewenvironment{cfg}{\lstset{language=cfg}}{}

% ------------------------------------------------------------------
% Lists
% Pyre inventory
\newcommand{\property}[2]{\item[\ttfamily\bfseries\color{purple}#1]{#2}}
\newcommand{\facility}[2]{\item[\ttfamily\bfseries\color{blue}#1]{#2}}
\newenvironment{inventory}{%
  \begin{list}{}{%
      \setlength{\topsep}{0em}%
      \setlength{\labelwidth}{1in}%
      \setlength{\leftmargin}{\labelwidth+\labelsep}%
      \setlength{\itemsep}{0.2em}%
      \setlength{\parsep}{0em}}}{%
  \end{list}
}

% End of file
