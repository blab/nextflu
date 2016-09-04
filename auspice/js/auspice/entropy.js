width = parseInt(d3.select(".entropy-container").style("width"), 10);
height = 250;
var position = "inset";

d3.json(path + file_prefix + "entropy.json", function(error, json){
	console.log('Loading entropy');
	console.log(error, path+file_prefix+"entropy.json");

	//gt_chart.axis.x.values = ticks;
	/**
		parses a genotype string into region and positions
	**/

	var chart_data = {}
	var chart_types = {}
	var chart_xaxis = {}
	var posToAA = {};
	var ymin = 0;
	var xmax = 0;
	var anno_count=0;
	entropy = json;
	console.log(entropy);
    for (x in entropy){
		if (x!='nuc'){
	        start = entropy[x]['pos'][0]; end = entropy[x]['pos'][entropy[x]['pos'].length-1];
	        chart_data['x'+x+'anno'] = [start, 0.5*(start+end), end];
	        chart_data[x+'anno'] = [anno_count%3, anno_count%3, anno_count%3].map(function (d) {return -0.1*(d+1);});
	        anno_count+=1;
	        if (ymin>chart_data[x+'anno'][0]){
	            ymin = chart_data[x+'anno'][0];
	        }
	        chart_types[x+'anno'] = 'line';
	        chart_xaxis[x+'anno'] = 'x'+x+'anno';
	    }
    }

	for (gene in entropy){
		if (gene!='nuc'){
			chart_data[gene]=[];
			chart_data['x'+gene]=[];
			chart_types[gene]='bar';
			chart_xaxis[gene]='x'+gene;
			for (var ii=0;ii<entropy[gene]["pos"].length;ii+=1){
				if (Math.round(10000*entropy[gene]["val"][ii])/10000>0.05){
					chart_data[gene].push(Math.round(10000*entropy[gene]["val"][ii])/10000);
					chart_data['x'+gene].push(entropy[gene]["pos"][ii]);
					posToAA[entropy[gene]["pos"][ii]] = [gene, entropy[gene]["codon"][ii]]
				}
			}
			var tmp_xmax = d3.max(chart_data['x'+gene]);
			if (tmp_xmax>xmax) {xmax=tmp_xmax;}
		}
	}
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
					values: ([1,2,3,4,5]).map(function (d){
						var dec = Math.pow(10,Math.floor(Math.log10(xmax/5)))
						var step = dec*Math.floor(xmax/5/dec);
						return d*step;
					})
				}
			},
		},
		data: {
			xs: chart_xaxis,
			json: chart_data,
			types: chart_types,
			onclick: function (d,i) {
            	gene = posToAA[d.x][0];
            	var pos = posToAA[d.x][1];
				if (entropy[gene]["val"][pos]>0.0){
					colorBy = "genotype";
					console.log("color by genotype: "+gene + ' ' + pos)
					colorByGenotypePosition([[gene, pos]]);
					d3.select("#gt-color").property("value", gene + ':' + (pos+1));
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
        	},
        	focus:{
        		show:false
        	}
    	},
	    tooltip: {
	        format: {
	            title: function (d) {
	            	if (typeof posToAA[d] != "undefined"){
		            	var gene = posToAA[d][0];
		            	var pos = posToAA[d][1];
		            	return gene + ' codon ' + (pos+1) +": " + entropy[gene]["val"][pos];
		            }else{ return d;}},
	            value: function (value, ratio, id) {
	                return id.substring(id.length-4)=='anno'?"start/stop":"Variability: "+value;
	            }
	        }
		},
	});
});
