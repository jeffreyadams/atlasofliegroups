\def\currentCWEBxversion{x3.51} % see also version_string end of this file

\def\hang{\hangindent 3em\indent\ignorespaces}
\def\pb{$\.|\ldots\.|$} % C brackets (|...|)
\def\LKC.{Levy/Knuth \.{CWEB}}
\def\:#1{`\.{@@#1}'}

\def\title{\me. (Version \currentCWEBxversion)}
\def\topofcontents
{\topglue 0pt plus .5 fill
 \centerline{\titlefont The {\ttitlefont \me.} program}
 \vskip 15pt
 \centerline{(\.{CWEB} version \currentCWEBxversion)}
 \vfill
}
\def\botofcontents
{\vfill\noindent
 Copyright \copyright\ 1987,\thinspace1990 Silvio Levy and Donald E. Knuth
 \par\noindent
 Copyright 1994 Marc A. A. van Leeuwen
 \bigskip\noindent
 Permission is granted to make and distribute verbatim copies of this
 document provided that the copyright notice and this permission notice
 are preserved on all copies.

 \smallskip\noindent
 Permission is granted to copy and distribute modified versions of this
 document under the conditions for verbatim copying, provided that the
 entire resulting derived work is distributed under the terms of a
 permission notice identical to this one.
}


@* Introduction. @:title@>
This is the main source text for the program \.{\me.}, one of the two text
processing tools that constitute the \.{CWEB} system. This version of the
system, identified as \.{CWEB}\currentCWEBxversion, was written by Marc van
Leeuwen. The history of this system is somewhat complicated, and starts with
the \.{WEB} system created by D.~E. Knuth (for the programming language
Pascal); the program corresponding to the current program in that system was
called~\myroots. That program was converted into \.{\me.}, written in and for
the \Cee~language, by Silvio Levy, and further developed under joint
responsibility with Knuth. Based on \.{CWEB}~2.1 by Levy and Knuth, a version
adapted to \caps{ANSI~C} was made by Frank Jensen. That version has served as
the starting point for the current author; an intermediate version between it
and the current version was made public under the name \.{CWEB}~3.x. Helpful
comments and suggestions by Gareth McCaughan, Werner Lemberg, and Andreas
Scherer were greatly appreciated.

The name \.{CWEB}~3.x was an unfortunate choice, since the Levy/Knuth
version of \.{CWEB} had also been independently developed further, and it is
currently being distributed as \.{CWEB}~3.4f. Both these branches have
retained the basic functionality of the version of \.{CWEB} they were both
derived from, making only minor (but different) additions and changes.
Realising that this divergence is undesirable and confusing, the development
of this branch from \.{CWEB}~3.x to the current version was mainly concerned
with improving compatibility and peaceful coexistence with the other branch.
First of all the name was changed by moving the character~`x' to the front
of the version number, so that this branch of the system can be referred to
as `\.{CWEBx}' without mentioning a complete version number; for clarity we
shall refer to the other branch as `\LKC.'. @^Levy/Knuth \.{CWEB}@>
Secondly, all the extensions of \LKC. were included, although in a few cases
in a slightly different form to avoid name conflicts with features already
present in \.{CWEB}~3.x. Thirdly, in order to accomodate existing programs
written for \LKC., among which are those of Knuth's impressive collection
called the Stanford GraphBase, a ``compatibility mode'' was added,
selectable by a command line option (\.{+c}), in which the remaining, mostly
trivial, differences in syntax and semantics are removed, and the system
attempts to be an alternative implementation for exactly the same language
as accepted by \LKC., at the price of losing a few possibilities particular
to~\.{CWEBx}.

The major version number `3' now corresponds to the same version number of
\LKC., and since no further extension of the possibilities of that system
are anticipated, we do not expect that compatibility will ever require any
further increase of this major version number. The Stanford GraphBase
@^Stanford GraphBase@> states that it requires \.{CWEB} version 3.0 or
greater; this version of \.{CWEBx} can fully process it in a satisfactory
way, using compatibility mode.

The ``banner line'' defined here should be changed whenever changes to
\.{\me.} are made publicly available. Users who are discontent with some
aspect of the program are encouraged to make the necessary modifications in
their copy of the source files; if this involves an improvement of the
implementation, they are kindly requested to inform the author responsible
for the current version. If they choose to make the resulting program
available under the name \.{\me.}, they should modify the non-numeric
prefix of the version specification.

@d version_string "x3.51"
