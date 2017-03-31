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
\RequirePackage{enumitem}
\RequirePackage[table]{xcolor}
\RequirePackage{listings}
\RequirePackage{pifont}
\RequirePackage{calc}
\RequirePackage{tikz}

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
\definecolor{mdblue}{rgb}{0.28, 0.50, 0.80} % 72/128/205
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
      ]{\raisebox{1em}{{\LARGE\bomb} {\bf Warning}}\\ #1};
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
\lstnewenvironment{cfg}{\lstset{language=cfg}}{}
%
% shell  
\lstnewenvironment{shell}{\lstset{%
    language=bash,%
    frame=single,%
    basicstyle=\small\ttfamily,%
    morekeywords={$$},%
    commentstyle=\itshape\color{orange},
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
  basicstyle=\small\ttfamily,
  frame=tb,%
  alsoletter={-},%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={mesh,dimension,use-index-zero,vertices,count,coordinates,cells,num-corners,simplices,material-ids,group,name,type,indices},
}
\lstnewenvironment{MeshIOAscii}{\lstset{language=MeshIOAscii}}{}
%
% SimpleIOAscii
\lstdefinelanguage{SimpleIOAscii}{%
  basicstyle=\small\ttfamily,
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
  basicstyle=\small\ttfamily,
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
  basicstyle=\small\ttfamily,
  frame=tb,%
  alsoletter={-},%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={TimeHistory,num-points,time-units},
}
\lstnewenvironment{TimeHistoryIO}{\lstset{language=TimeHistoryIO}}{}
%
% PointsList
\lstdefinelanguage{PointsList}{%
  basicstyle=\small\ttfamily,
  frame=tb,%
  morecomment=[l][\color{orange}\itshape]{\#},
  morekeywords={units},
}
\lstnewenvironment{PointsList}{\lstset{language=PointsList}}{}
%
% TimeStepUser
\lstdefinelanguage{TimeStepUser}{%
  basicstyle=\small\ttfamily,
  frame=tb,%
  morecomment=[l][\color{orange}\itshape]{//},
  morekeywords={units},
}
\lstnewenvironment{TimeStepUser}{\lstset{language=TimeStepUser}}{}
%
% % Features table
\newcommand{\rlabel}[1]{\rotatebox[origin=l]{90}{#1}}
\newcommand{\yes}{\ding{52}}


% End of file
