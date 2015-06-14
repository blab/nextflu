var frequencies, pivots;
var gene = 'nuc'
/**
 * for each node, calculate the derivative of the frequency tranjectory. if none exists, copy parent
**/
function calcDfreq(node, freq_ii){
	if (typeof node.children != "undefined") {
		for (var i1=0; i1<node.children.length; i1++) {
			if (typeof node.children[i1].freq != "undefined") {
				if (node.children[i1].freq["global"] != "undefined"){
					var tmp_freq = node.children[i1].freq["global"]
					node.children[i1].dfreq = 0.5*(tmp_freq[freq_ii] - tmp_freq[freq_ii-dfreq_dn])/(tmp_freq[freq_ii] + tmp_freq[freq_ii-dfreq_dn] + 0.1);
				} else {
					node.children[i1].dfreq = node.dfreq;
				}
			}
			calcDfreq(node.children[i1], freq_ii);
		}
	}
};
function parse_gt_string(gt){
	mutations = [];
	gt.split(',').map( function (d) {
		var tmp = d.split(/[\s//]/); //FIXME: make more inclusive
		var region;
		var positions = [];
		for (var i=0; i<tmp.length; i++){
			if (contains(["EU","NA","AS","OC"], tmp[i])){
				region = tmp[i];
			}else{
				if (tmp[i].length>0) positions.push(tmp[i]);
			}
		}
		if (typeof region == "undefined") region="global";
		// sort of this is a multi mutation genotype
		if (positions.length>1){
			positions.sort(function (a,b){
				return parseInt(a.substring(0,a.length-1)) - parseInt(b.substring(0,b.length-1));
			});
		}
		mutations.push([region, positions.join('/')]);
	});
	return mutations;
};

/**
loops over all genotypes from a certain region and sums the frequency contributions
of the genotype matches at the specified positions
**/
function get_frequencies(region, gt){
	var freq = [];
	for (var pi=0; pi<pivots.length; pi++){freq[freq.length]=0;}
	if (frequencies["clades"][region][gt.toLowerCase()]!=undefined) {
		console.log(gt+" found as clade");
		for (var pi=0; pi<freq.length; pi++){
			freq[pi]+=frequencies["clades"][region][gt.toLowerCase()][pi];
		}
	}
	else if ((typeof frequencies["genotypes"] !="undefined") && (frequencies["genotypes"][region][gt]!=undefined)) {
		console.log(gt+" found as genotype");
		for (var pi=0; pi<freq.length; pi++){
			freq[pi]+=frequencies["genotypes"][region][gt][pi];
		}
	}else if (frequencies["mutations"][region][gt]!=undefined) {
		console.log(gt+" found as mutation");
		for (var pi=0; pi<freq.length; pi++){
			freq[pi]+=frequencies["mutations"][region][gt][pi];
		}
	}
	return freq.map(function (d) {return Math.round(d*100)/100;});
};


function make_gt_chart(gt){
	var tmp_data = [];
	var tmp_trace = ['x'];
	var tmp_colors = {};
	tmp_data.push(tmp_trace.concat(pivots));
	gt.forEach(function (d, i) {
		var region = d[0];
		var genotype = d[1];
		var freq = get_frequencies(region, genotype);
		if (d3.max(freq)>0) {
			var tmp_trace = genotype.toString().replace(/,/g, ', ');
			if (region != "global") {
				tmp_trace = region + ':\t' + tmp_trace;
			}
			tmp_data.push([tmp_trace].concat(freq));
			tmp_colors[tmp_trace] = genotypeColors[i];
		}
	});
	console.log(tmp_colors);
	gt_chart.load({
       	columns: tmp_data,
       	unload: true
	});
	gt_chart.data.colors(tmp_colors);
}

function addClade(d) {
	if (typeof gt_chart != "undefined"){
		var plot_data = [['x'].concat(rootNode["pivots"])];
		var reg = "global";
		if ((typeof d.target.freq !="undefined" )&&(d.target.freq[reg] != "undefined")){
			plot_data[plot_data.length] = [reg].concat(d.target.freq[reg]);				
		}
		if (plot_data.length > 1) {
			if (plot_data[1][0] == "global") {
				plot_data[1][0] = "clade";
			}
		}
		gt_chart.load({
	       	columns: plot_data
		});
	}
}

function removeClade() {
	if (typeof gt_chart != "undefined"){
		gt_chart.unload({
	       	ids: ["clade"]
		});
	}
}

width = parseInt(d3.select(".freqplot-container").style("width"), 10);
height = 250;
var position = "inset";

var gt_chart = c3.generate({
	bindto: '#gtchart',
	size: {width: width-10, height: height},
	onresize: function() {
		width = parseInt(d3.select(".freqplot-container").style("width"), 10);
		height = 250;
		gt_chart.resize({height: height, width: width});
	},
	legend: {
		position: position,
		inset: {
    		anchor: 'top-right',
    		x: 10,
    		y: -15,
    		step: 1
    	}
	},
  	color: {
        pattern: genotypeColors
    },
	axis: {
		y: {
			label: {
				text: 'frequency',
				position: 'outer-middle'	
			},
			tick: {
				values: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
				outer: false
			},
            min: 0,			
			max: 1	
		},
		x: {
			label: {
				text: 'time',
				position: 'outer-center'	
			},
			tick: {
				values: time_ticks,
				outer: false				
			}
		}
	},			
	data: {
		x: 'x',
		columns: [],
	}
});

function contains(arr, obj) {
    for(var i=0; i<arr.length; i++) {
        if (arr[i] == obj) return true;
    }
}

var case_numbers;
d3.json(path + file_prefix + "case_numbers.json", function(error, json){
	console.log(error);
	case_numbers = json;
	var case_numbers_chart = c3.generate({
		bindto: '#casenumberschart',
		size: {width: width-10, height: height},
		onresize: function() {
			width = parseInt(d3.select(".freqplot-container").style("width"), 10);
			height = 250;
			case_numbers_chart.resize({height: height, width: width});
		},
		legend: {
			position: position,
			inset: {
	    		anchor: 'top-right',
	    		x: 10,
	    		y: -15,
	    		step: 4
	    	}
		},
		axis: {
			y: {
				label: {
					text: 'weekly case numbers',
					position: 'outer-middle'	
				},
				tick:{
					values: [0, 50, 100, 150, 200, 250]
				},
			},
			x: {
				label: {
					text: 'time',
					position: 'outer-center'	
				},
				tick: {
					values: time_ticks,
					outer: false				
				}
			}
		},			
		bar: {width:3},
		data: {
			json: case_numbers,
			type: 'line',
			xs: {
				'Korea/China':'xKorea/China',
				'other':'xother',
				'Saudi Arabia':'xSaudi Arabia'
			}
		}
	});
});


d3.json(path + file_prefix + "frequencies.json", function(error, json){
	console.log(error);
	frequencies = json;
	pivots= frequencies["mutations"]["global"]["pivots"].map(function (d) {return Math.round(parseFloat(d)*100)/100;});
	var ticks = [Math.round(pivots[0])];
	var tick_step = Math.round((pivots[pivots.length-1]-pivots[0])/6*10)/10;
	while (ticks[ticks.length-1]<pivots[pivots.length-1]){
		ticks.push(Math.round((ticks[ticks.length-1]+tick_step)*10)/10);
	}
	//gt_chart.axis.x.values = ticks;
	/**
		parses a genotype string into region and positions
	**/

	var chart_data = {}
	var chart_types = {}
	var chart_xaxis = {}
	var posToAA = {};
	var ymin = 0;
	if (typeof genome_annotation !== 'undefined') {
		for (x in genome_annotation){
			chart_data['x'+x+'anno'] = genome_annotation[x][1];
			chart_data[x+'anno'] = genome_annotation[x][0].map(function(d) {return -0.1*d;});
			if (ymin>chart_data[x+'anno'][0]){
				ymin = chart_data[x+'anno'][0];
			}
			chart_types[x+'anno'] = 'line';
			chart_xaxis[x+'anno'] = 'x'+x+'anno';
		}
		ymin-=0.08;
	}

	for (gene in frequencies["entropy"]){
		console.log('adding gene to entropy: '+gene);
		chart_data[gene]=[];
		chart_data['x'+gene]=[];
		chart_types[gene]='bar';
		chart_xaxis[gene]='x'+gene;
		if 	(typeof chart_data['x'+gene+'anno']!="undefined"){
			var offset= chart_data['x'+gene+'anno'][0];
		}else{
			var offset=1;
		}
		for (var ii=0;ii<frequencies["entropy"][gene].length;ii+=1){
			if (Math.round(10000*frequencies["entropy"][gene][ii][1])/10000>0.05){
				chart_data[gene].push(Math.round(10000*frequencies["entropy"][gene][ii][1])/10000);
				chart_data['x'+gene].push(ii*3+offset);
				posToAA[ii*3+offset] = [gene, ii];
			}
		}
	}
	console.log(chart_data);
	console.log(chart_types);
	console.log(chart_xaxis);
	var entropy_chart = c3.generate({
		bindto: '#entropy',
		size: {width: width-10, height: height},
		onresize: function() {
			width = parseInt(d3.select(".entropy-container").style("width"), 10);
			height = 250;
			entropy_chart.resize({height: height, width: width});
		},		
		legend: {show: false},
		color: {pattern: ["#AAA"]},
		axis: {
			y: {
				label: {
					text: 'variability',
					position: 'outer-middle'	
				},
				tick: {
					values: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
					outer: false
				},
				min:ymin,
			},
			x: {
				label: {
					text: 'position',
					position: 'outer-center'	
				},
				tick: {
					outer: false,
					values:[0, 5000, 10000, 15000, 20000, 25000, 30000]
				}
			},
		},			
		data: {
			xs: chart_xaxis,
			json: chart_data,
			types: chart_types,
			onclick: function (d,i) { 
            	var gene = posToAA[d.x][0];
            	var pos = posToAA[d.x][1];
				console.log(d);
				if (frequencies["entropy"][gene][pos][2].length>1){
					var tmp = [];
					for (var ii=0;ii<frequencies["entropy"][gene][pos][2].length;ii+=1){
						tmp.push(["global",d.x+frequencies["entropy"][gene][pos][2][ii]]);
					}
					colorBy = "genotype";
					colorByGenotypePosition([pos], gene);
					d3.select("#gt-color").property("value", gene + ' ' + pos);
				}
		    },
		    onmouseover: function (d){
		    	document.body.style.cursor = "pointer";
		    },
		    onmouseout: function (d){
		    	document.body.style.cursor = "default";
		    },
			labels:{
				format:function (v, id, i, j){
					if ((typeof id !="undefined")&&(id.substring(id.length-4)=='anno')&&(i==1)){
						return id.substring(0,id.length-4);
					}else{return '';}
				}
			},
		},
		bar: {width: 2},
	    grid: {
    	    y: {
        	    lines: [{value: 0}]
        	}
    	},
	    tooltip: {
	        format: {
	            title: function (d) {
	            	if (typeof posToAA[d] != "undefined"){
		            	var gene = posToAA[d][0];
		            	var pos = posToAA[d][1];
		            	console.log(gene+ ' '+pos);
		                {
							if (frequencies["entropy"][gene][pos][2].length>1){
								var tmp = [];
								for (var ii=0;ii<frequencies["entropy"][gene][pos][2].length;ii+=1){
									tmp.push(["global",d+frequencies["entropy"][gene][pos][2][ii]]);
								}
								console.log("tooltip");
								colorBy = "genotype";
								colorByGenotypePosition([pos], gene);
								d3.select("#gt-color").property("value", d);
							}
						}
		            	return 'Genomic position ' + d + ', ' + gene + ' codon ' + pos + frequencies["entropy"][gene][pos][2].join(",");
		            }else{ return d;}},
	            value: function (value, ratio, id) {
	                return id.substring(id.length-4)=='anno'?"start/stop":"Variability: "+value;
	            }
	        }
		},
	});

	d3.select("#plotfreq")
		.on("click", function (){
			gt = parse_gt_string(document.getElementById("gtspec").value);			
			make_gt_chart(gt);
		});
	make_gt_chart(parse_gt_string(document.getElementById("gtspec").value));
});
