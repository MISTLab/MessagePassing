# MessagePassing

## Requirements
Folder `buzz/` requires
```
argos3   # see https://www.argos-sim.info/
buzzc   # see https://github.com/MISTLab/Buzz
```

Folder `spectral-graph-theory/` requires
```
octave   # https://www.gnu.org/software/octave/
```

## How to use
Folder `buzz/`
```
git clone https://github.com/MISTLab/MessagePassing
cd MessagePassing/buzz/
bzzc control.bzz
argos3 -c sim.argos
```

Folder `spectral-graph-theory/`
```
git clone https://github.com/MISTLab/MessagePassing
cd MessagePassing/spectral-graph-theory/
octave runAll.m
cd latex
pdflatex plot.tex
```
