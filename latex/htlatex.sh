#!/bin/bash

# Compile the LaTeX document using xelatex
xelatex report.tex

# Run htlatex if needed
htlatex report.tex

