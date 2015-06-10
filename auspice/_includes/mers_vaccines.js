var	vaccineChoice = {};
var vaccineStrains = Object.keys(vaccineChoice);

var genome_annotation = {'ORF1a':[[1,1,1], [233, 6500, 13408]],
						 'ORF1b':[[2,2,2], [13387, 15000,  21468]],
						 'S':[[1,1,1], [21410, 23000, 25471]],
						 'ORF3':[[3,3,3], [ 25486, 25500, 25797]],
						 'ORF4ab':[[1,1,1], [ 25806, 25900, 26787]],
						 'ORF5':[[3,3,3], [ 26794, 27000, 27468]],
						 'E':[[1,1,1], [ 27544, 27600, 27792]],
						 'M':[[3,3,3], [ 27807, 28000, 28466]],
						 'N':[[1,1,1], [ 28520, 29000, 29761]],
						 'ORF8b':[[3,3,3], [ 28716, 29054, 29054]],
}
var restrictTo = {"country":"all","region":"all", "lab":"all", "host":"all"};

var codeToCountry = {
"ARE":"U. Arab Emirates",
"EGY":"Egypt",
"KSA":"Saudi Arabia",
"JOR":"Jordan",
"GBR":"United Kingdom",
"QAT":"Qatar",
"OMN":"Oman",
"China":"China",
"KOR":"South Korea"
}