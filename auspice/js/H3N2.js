---
---
var file_prefix = 'H3N2_';
var dfreq_dn = 2;

var time_ticks = [2012.5, 2013, 2013.5, 2014, 2014.5, 2015];
var	time_window = 1.0;  // layer of one year that is considered current or active
var LBItau = 0.0008;
var LBItime_window = 0.5;
var freqdefault = "3c2.a, 3c3.a"

{%include_relative H3N2_vaccines.js %}
