% tgammazcrossing.mp -- do not edit, generated automatically by Report.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 40*2.84526pt#, 25*2.84526pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=2.84526pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1, __i2);
vright (__o1, __o2);
vlabel ("$e^+, p_{2}$", __i1);
vlabel ("$e^-, p_{1}$", __i2);
vlabel ("$\mu ^+, p_{4}$", __o1);
vlabel ("$\mu ^-, p_{3}$", __o2);
vconnect ("fermion", __o1, __v1, __i1);
vconnect ("fermion", __i2, __v2, __o2);
vconnect ("photon,label=$\gamma ^{*}/\Zzero ^{(*)}$", __v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
