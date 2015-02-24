## Introduction

nextflu is designed to perform near real-time tracking of influenza virus evolution. It's divided into two components: [augur](augur/), which takes a `.fasta` file of flu sequences and builds an annotated phylogeny, and [auspice](auspice/), which displays this annotated phylogeny in an interactive web-based visualization.

The current version of nextflu is focused on tracking seasonal influenza H3N2 evolution in humans, looking at sequences from the hemagglutinin (HA) gene. Future versions may extend this analysis to other genes in H3N2 or other influenza subtypes. We would also like to implement formal predictive models to make nextflu a platform for forecasting evolution in addition to up-to-date tracking.

## License

All source code in this repository is freely available under an MIT license, unless otherwise noted within a file.

**The MIT License (MIT)**

Copyright (c) 2015 Trevor Bedford and Richard Neher

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.