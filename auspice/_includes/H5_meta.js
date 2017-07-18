var	vaccineChoice = {};
var vaccineStrains = Object.keys(vaccineChoice);
var subtypes = ['Unknown', 'H5N1','H5N2','H5N3','H5N4', 'H5N5','H5N6','H5N7','H5N8','H5N9'];



var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [1,20,49]],
                         'HA1':[[1,1,1], [49,460,49+987]],
						 'HA2':[[1.2,1.2,1.2], [49+987,1200,1701]]}
