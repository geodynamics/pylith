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
% See LICENSE.md for license information.
%
% ----------------------------------------------------------------------
%
% pylithdoc.cls    Style file for PyLith documentation.
\NeedsTeXFormat{LaTeX2e}[1994/06/01]
\ProvidesClass{pylithdoc}[2017/09/01 v0.2 PyLith Documentation class]

\ProcessOptions*
\LoadClass[10pt,twoside]{book}

% Packages
\RequirePackage{fourier} % texlive-fonts-extra
\RequirePackage{float}
\RequirePackage{varioref} % \vref
\RequirePackage{enumitem}
\RequirePackage[table]{xcolor}
\RequirePackage{listings}
\RequirePackage{pifont}
\RequirePackage{calc}
%\RequirePackage[table]{xcolor}
\RequirePackage{tikz}
\RequirePackage{makecell} % tables
\RequirePackage{booktabs} % tables
\RequirePackage{colortbl} % tables


%\RequirePackage{url}
\RequirePackage[%
  unicode=true,%
  pdfusetitle,%
  pdftex,%
  bookmarks=true,%
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
\RequirePackage{helvet}
\renewcommand{\LARGE}{\@setfontsize\LARGE{17}{22}}
\renewcommand{\Large}{\@setfontsize\Large{14.4}{18}}
\renewcommand{\large}{\@setfontsize\large{10.95}{13.6}}
\renewcommand{\normalsize}{\@setfontsize\normalsize{10}{12}}
\renewcommand{\small}{\@setfontsize\small{9}{10}}
\renewcommand{\footnotesize}{\@setfontsize\footnotesize{8}{9.5}}
%
% Colors
\definecolor{ltred}{rgb}{1.0, 0.57, 0.59} % 255/145/151
\definecolor{red}{rgb}{0.63, 0.00, 0.06} % 161/0/16

\definecolor{ltblue}{rgb}{0.57, 0.59, 1.00} % 145/151/255
\definecolor{mdblue}{rgb}{0.28, 0.50, 0.80} % 72/128/205
\definecolor{blue}{rgb}{0.00, 0.06, 0.63} % 0/16/161

\definecolor{ltgreen}{rgb}{0.59, 1.00, 0.57} % 151/255/145
\definecolor{green}{rgb}{0.06, 0.63, 0.00} % 16/161/0

\definecolor{ltpurple}{rgb}{0.62, 0.29, 0.63} % 157/75/161
\definecolor{purple}{rgb}{0.63, 0.00, 0.54} % 131/1/138

\definecolor{ltcyan}{rgb}{0.33, 0.70, 0.68} % 83/179/174
\definecolor{cyan}{rgb}{0.00, 0.54, 0.51} % 1/138/131

\definecolor{ltyellow}{rgb}{1.0, 0.97, 0.46} % 255/248/118
\definecolor{yellow}{rgb}{0.90, 0.83, 0.0} % % 230/211/0


\definecolor{ltorange}{rgb}{1.0, 0.74, 0.35} % 255/189/89
\definecolor{orange}{rgb}{0.87, 0.45, 0.0} % 222/115/0

\definecolor{dkslate}{rgb}{0.18, 0.21, 0.28} % 47/53/72
\definecolor{mdslate}{rgb}{0.45, 0.50, 0.68} % 114/127/173
\definecolor{ltslate}{rgb}{0.85, 0.88, 0.95} % 216/225/229

% ------------------------------------------------------------------
% Headings
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{3}

% ------------------------------------------------------------------
% Figures
% automatic centering
\g@addto@macro\@floatboxreset\centering

% ------------------------------------------------------------------
% Macros
%
% PyLith version
\newcommand{\pylithVersionNumber}{X.X.X}
\newcommand{\pylithDOI}{10.5281/zenodo.XXXXXX}
\newcommand{\pylithVersion}{v\pylithVersionNumber}
% Boxes
\usetikzlibrary{shadows.blur}
% Warning
\newcommand{\userwarning}[1]{%
  \begin{center}\begin{tikzpicture}
      \node[draw=none,
        shade,
        text width=\textwidth-2.0in,
        inner sep=1em,
        top color=red!50,
        bottom color=red!5,
        rounded corners=1em,
        blur shadow={shadow blur steps=5}
      ]{\raisebox{1em}{{\LARGE\bomb} {\bf Warning}}\\ #1};
    \end{tikzpicture}\end{center}}
%
% Tip
\newcommand{\usertip}[1]{%
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
      ]{\raisebox{1em}{{\LARGE\danger} {\bf Important}}\\ #1};
    \end{tikzpicture}\end{center}}
%
% TODO
\newcommand{\todo}[2]{%
  \begin{center}\begin{tikzpicture}
      \node[draw=none,
        shade,
        text width=\textwidth-2.0in,
        inner sep=1em,
        top color=yellow,
        bottom color=yellow!10,
        rounded corners=1em,
        blur shadow={shadow blur steps=5}
      ]{\raisebox{1em}{{\LARGE\noway} {\bf TODO @#1}}\\ #2};
    \end{tikzpicture}\end{center}}
\newcommand{\todocomment}[2]{{\bfseries @#1} #2}
%
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
% file, directory
\newcommand{\filename}[1]{{\ttfamily #1}}
\newcommand{\commandline}[1]{{\ttfamily #1}}
\newcommand{\object}[1]{{\sffamily #1}}
\newcommand{\menu}[1]{{\sffamily #1}}
\newcommand{\gitbranch}[1]{{\ttfamily #1}}

% ------------------------------------------------------------------
% Code blocks
%
% Caption
\renewcommand{\lst@makecaption}[2]{%
  \if\lst@title\@empty
  \else
  \vskip\abovecaptionskip%
  \begin{tikzpicture}
    \node[draw=none,
      inner sep=0.3em,
      top color=\captioncolor!50,
      bottom color=\captioncolor!10,
    ]{#1: #2};
  \end{tikzpicture}\par\vspace{-1pt}
  \fi
}
%
%
% .cfg
%
\lstdefinelanguage{cfg}{%
  basicstyle=\footnotesize\ttfamily,
  frame=tb,%
  morecomment=[l][\color{orange}\itshape]{\#},
  morecomment=[l][\color{orange}\itshape]{;},
  moredelim=[is][keywordstyle]{<h>}{</h>},
  moredelim=[is][\color{purple}\bfseries]{<p>}{</p>},
  moredelim=[is][\color{blue}\bfseries]{<f>}{</f>},
}
\lstnewenvironment{cfg}[1][\@empty]{%
  \def\captioncolor{cyan}
  \def\@makecaption{\lst@makecaption}
  \lstset{%
    language=cfg,
    title={#1},
    backgroundcolor=\color{\captioncolor!10},
    rulecolor=\color{\captioncolor},
}}{}
%
% shell
%
\lstnewenvironment{shell}[1][\@empty]{%
  \def\captioncolor{mdslate}
  \def\@makecaption{\lst@makecaption}
  \lstset{%
    language=bash,%
    frame=tb,%
    basicstyle=\footnotesize\ttfamily,%
    morekeywords={$},%
    commentstyle=\itshape\color{orange},
    title={#1},
    backgroundcolor=\color{\captioncolor!10},
    rulecolor=\color{\captioncolor},
}}{}
%
% python
%
\lstnewenvironment{python}[1][\@empty]{%
  \def\captioncolor{purple}
  \def\@makecaption{\lst@makecaption}
  \lstset{%
    language=python,%
    frame=tb,%
    basicstyle=\footnotesize\ttfamily,%
    morekeywords={>>>},%
    commentstyle=\itshape\color{orange},
    title={#1},
    backgroundcolor=\color{\captioncolor!10},
    rulecolor=\color{\captioncolor},
}}{}
%
% C++
%
\lstnewenvironment{cplusplus}[1][\@empty]{%
  \def\captioncolor{ltgreen}
  \def\@makecaption{\lst@makecaption}
  \lstset{%
    language=c++,%
    frame=tb,%
    basicstyle=\footnotesize\ttfamily,%
    commentstyle=\itshape\color{orange},
    title={#1},
    backgroundcolor=\color{\captioncolor!10},
    rulecolor=\color{\captioncolor},
}}{}

% ------------------------------------------------------------------
% Lists
%
\setlist[itemize]{itemsep=0.1em,parsep=0em}
% Pyre inventory
\newcommand{\propertyshape}[1]{{\ttfamily #1}}
\newcommand{\facilityshape}[1]{{\ttfamily #1}}
\newcommand{\property}[1]{{\ttfamily\bfseries\color{purple}#1}}
\newcommand{\facility}[1]{{\ttfamily\bfseries\color{blue}#1}}
\newcommand{\propertyitem}[2]{\item[\property{#1}]{#2}}
\newcommand{\facilityitem}[2]{\item[\facility{#1}]{#2}}
\newenvironment{inventory}{%
  \begin{list}{}{%
      \setlength{\topsep}{0em}%
      \setlength{\labelwidth}{1in}%
      \setlength{\leftmargin}{\labelwidth+\labelsep}%
      \setlength{\itemsep}{0.2em}%
      \setlength{\parsep}{0em}}}{%
  \end{list}
}

% ------------------------------------------------------------------
% File formats
% MeshIOAscii
\lstdefinelanguage{MeshIOAscii}{%
  basicstyle=\footnotesize\ttfamily,
  frame=tb,%
  alsoletter={-},%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={mesh,dimension,use-index-zero,vertices,count,coordinates,cells,num-corners,simplices,material-ids,group,name,type,indices},
}
\lstnewenvironment{MeshIOAscii}{\lstset{language=MeshIOAscii}}{}
%
% SimpleIOAscii
\lstdefinelanguage{SimpleIOAscii}{%
  basicstyle=\footnotesize\ttfamily,
  frame=tb,%
  alsoletter={-},%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={SimpleDB,num-values,value-units,value-names,space-dim,cs-data,num-locs,data-dim},
  morekeywords={to-meters,ellipsoid,datum-horiz,datum-vert,is-geocentric,origin-lon,origin-lat,origin-elev,rotation-angle,projector,projection,units,proj-options},
}
\lstnewenvironment{SimpleIOAscii}{\lstset{language=SimpleIOAscii}}{}
%
% SimpleGridDB
\lstdefinelanguage{SimpleGridDB}{%
  basicstyle=\footnotesize\ttfamily,
  frame=tb,%
  alsoletter={-},%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={SimpleGridDB,num-values,value-units,value-names,space-dim,cs-data,num-locs,num-x,num-y,num-z,data-dim},
  morekeywords={to-meters,ellipsoid,datum-horiz,datum-vert,is-geocentric,origin-lon,origin-lat,origin-elev,rotation-angle,projector,projection,units,proj-options},
}
\lstnewenvironment{SimpleGridDB}{\lstset{language=SimpleGridDB}}{}
%
% TimeHistoryIO
\lstdefinelanguage{TimeHistoryIO}{%
  basicstyle=\footnotesize\ttfamily,
  frame=tb,%
  alsoletter={-},%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={TimeHistory,num-points,time-units},
}
\lstnewenvironment{TimeHistoryIO}{\lstset{language=TimeHistoryIO}}{}
%
% PointsList
\lstdefinelanguage{PointsList}{%
  basicstyle=\footnotesize\ttfamily,
  frame=tb,%
  morecomment=[l][\color{orange}\itshape]{\#},
  morekeywords={units},
}
\lstnewenvironment{PointsList}{\lstset{language=PointsList}}{}
%
% TimeStepUser
\lstdefinelanguage{TimeStepUser}{%
  basicstyle=\footnotesize\ttfamily,
  frame=tb,%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={units},
}
\lstnewenvironment{TimeStepUser}{\lstset{language=TimeStepUser}}{}
%
% % Features table
\newcommand{\rlabel}[1]{\rotatebox[origin=l]{90}{#1}}
\newcommand{\yes}{\ding{52}}
%
% Tables
\renewcommand{\thead}[1]{{\bfseries #1}}

% End of file
