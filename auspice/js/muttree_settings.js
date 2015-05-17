var file_prefix = '';
var path = '';
var dfreq_dn = 2;
var branch_labels = true;
var tip_labels = true;

var	vaccineChoice = {};
var vaccineStrains = Object.keys(vaccineChoice);
var genericDomain = [0,0.111,0.222,0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 1.0];
var epiColorDomain = [4,5,6,7,8,9,10,11,12,13];
var nonEpiColorDomain = [2,3,4,5,6,7,8,9,10,11];
var dateColorDomain = genericDomain;
var rbsColorDomain = [0,1,2,3,4];
var dfreqColorDomain = genericDomain.map(function(d){return Math.round(100*(-0.18+d*0.36))/100;});

var time_ticks = [2012.5, 2013, 2013.5, 2014, 2014.5, 2015];
var	time_window;  // layer of one year that is considered current or active
var LBItau = 0.0008;
var LBItime_window = 0.5;

