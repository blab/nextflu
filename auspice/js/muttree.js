var file_prefix = '';

function gatherTips(node, tips) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			gatherTips(node.children[i], tips);
		}
	}
	else {
		tips.push(node);
	}
	return tips;
}

function gatherInternals(node, internals) {
	if (typeof node.children != "undefined") {
		internals.push(node);
		for (var i=0, c=node.children.length; i<c; i++) {
			gatherInternals(node.children[i], internals);
		}
	}
	return internals;
}

function minimumAttribute(node, attr, min) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			min = minimumAttribute(node.children[i], attr, min);
		}
	}
	else {
		if (node[attr] < min) {
			min = node[attr];
		}
	}
	return min;
}

function maximumAttribute(node, attr, max) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			max = maximumAttribute(node.children[i], attr, max);
		}
	}
	else {
		if (node[attr] > max) {
			max = node[attr];
		}
	}
	return max;
}

function tipRadius(d) {
	var radius = 4;
	return radius;
}

function branchStrokeColor(col) {
	var modCol = d3.interpolateRgb(col, "#BBB")(0.6);
	return d3.rgb(modCol).toString();
}

function tipStrokeColor(col) {
	return d3.rgb(col).toString();	
}

function tipFillColor(col) {
	return d3.rgb(col).brighter([0.65]).toString();
}

function treePlotHeight(width) {
	return 400 + 0.35*width;
}

var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);

var width = containerWidth,
	height = treePlotHeight(containerWidth);

var cladeToSeq = {}

var globalDate = new Date();
var ymd_format = d3.time.format("%Y-%m-%d");

var tree = d3.layout.tree()
	.size([height, width]);

var treeplot = d3.select("#treeplot")
	.attr("width", width)
	.attr("height", height);

var legend = d3.select("#legend")
	.attr("width", 280)
	.attr("height", 100);

var virusTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string='';
		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";
		
		string += "<div class=\"smallnote\">";		
		
		if (typeof d.country != "undefined") {
			string += d.country.replace(/([A-Z])/g, ' $1');
		}
		if (typeof d.date != "undefined") {
			string += ", " + d.date;
		}
		if (typeof d.isolate_id != "undefined") {
			string += "<br>GISAID ID: EPI" + d.accession;
		}
		if (typeof d.orig_lab != "undefined") {
			if (d.orig_lab != "") {
				string += "<br>Source: " + d.orig_lab.substring(0,25);
				if (d.orig_lab.length>25) string += '...';
			}
		}			
		if (typeof d.sub_lab != "undefined") {
			if (d.sub_lab != "") {
				string += "<br>Subm: " + d.sub_lab.substring(0,25);
				if (d.sub_lab.length>25) string += '...';
			}
		}			
		string += "</div>";
		console.log(d.desc);
		return string;
	});
treeplot.call(virusTooltip);

var linkTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = "Mutations: ";
		if (d.aa_muts.length){
			string+=d.aa_muts.replace(/,/g, ', ');
		}
		return string;
	});
treeplot.call(linkTooltip);
var colorBy = "genotype";
var colorScale;

d3.json(file_prefix + "tree.json", function(error, root) {

	if (error) return console.warn(error);

	var nodes = tree.nodes(root),
		links = tree.links(nodes);
	var tree_legend;
	var rootNode = nodes[0];
	var tips = gatherTips(rootNode, []);
	var internals = gatherInternals(rootNode, []);

	var xValues = nodes.map(function(d) {
		return +d.xvalue;
	});

	var yValues = nodes.map(function(d) {
		return +d.yvalue;
	});

	var xScale = d3.scale.linear()
		.domain([d3.min(xValues), d3.max(xValues)])
		.range([10, width-50]);

	var yScale = d3.scale.linear()
		.domain([d3.min(yValues), d3.max(yValues)])
		.range([20, height-10]);

	nodes.forEach(function (d) {
		d.x = xScale(d.xvalue);
		d.y = yScale(d.yvalue);
	});

	var recencyLinksSizeScale = d3.scale.threshold()
		.domain([0.0])
		.range([0, 2]);

	var genotypeColors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"]

    var legendRectSize = 15;
    var legendSpacing = 4;
    function makeLegend(){
    	
    	d3.select("#legend-title").text(function(d){
   			if (colorBy == "genotype") {
    			return "Genotype"
    		}
    	});
    
		var tmp_leg = legend.selectAll(".legend")
		    .data(colorScale.domain())
		    .enter().append('g')
		    .attr('class', 'legend')
		    .attr('transform', function(d, i) {
		    	var stack = 5;
				var height = legendRectSize + legendSpacing;
				var fromRight = Math.floor(i / stack);
				var fromTop = i % stack;
				var horz = fromRight * 145 + 5;				
				var vert = fromTop * height + 5;
				return 'translate(' + horz + ',' + vert + ')';
		    	});
		tmp_leg.append('rect')
		    .attr('width', legendRectSize)
		    .attr('height', legendRectSize)
		    .style('fill', function (d) {
		   		var col = colorScale(d);
		   		return d3.rgb(col).brighter([0.35]).toString();
		    })
		    .style('stroke', function (d) {
		    	var col = colorScale(d);
		    	return tipStrokeColor(col);
		    });
		
		tmp_leg.append('text')
		    .attr('x', legendRectSize + legendSpacing + 5)
		    .attr('y', legendRectSize - legendSpacing)
		    .text(function(d) {
		    	return d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');
		    });		
		return tmp_leg;
    }

    function removeLegend(){
    	legend.selectAll('.legend')
    		.remove();
    }


	var link = treeplot.selectAll(".link")
		.data(links)
		.enter().append("polyline")
		.attr("class", "link")
		.attr("points", function(d) {
			var mod = 1;
			return (d.source.x-mod).toString() + "," + d.source.y.toString() + " "
			+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
			+ (d.target.x).toString() + "," + d.target.y.toString()
		})
		.style("stroke-width", 3)
		.style("stroke", function(d) {
				var col = '#BBB';
				return branchStrokeColor(col);
			})		
		.style("cursor", "pointer")
		.on('mouseover', function(d) {
			linkTooltip.show(d.target, this);
		})
		.on('mouseout', linkTooltip.hide)		
		.on('click', function(d) {
			var dMin = 0.5 * (minimumAttribute(d.target, "xvalue", d.target.xvalue) + minimumAttribute(d.source, "xvalue", d.source.xvalue)),
				dMax = maximumAttribute(d.target, "xvalue", d.target.xvalue),
				lMin = minimumAttribute(d.target, "yvalue", d.target.yvalue),
				lMax = maximumAttribute(d.target, "yvalue", d.target.yvalue);
			if (dMax > dMin && lMax > lMin) {
				rescale(dMin, dMax, lMin, lMax);
			}
			else {
				dMin = minimumAttribute(d.source, "xvalue", d.source.xvalue),
				dMax = maximumAttribute(d.source, "xvalue", d.source.xvalue),
				lMin = minimumAttribute(d.source, "yvalue", d.source.yvalue),
				lMax = maximumAttribute(d.source, "yvalue", d.source.yvalue);			
				rescale(dMin, dMax, lMin, lMax);
			}
		});

	var tipCircles = treeplot.selectAll(".tip")
		.data(tips)
		.enter()
		.append("circle")
		.attr("class", "tip")
//		.attr("id", function(d) { return (d.strain).replace(/\//g, ""); })
		.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; })
		.attr("r", function(d) { return tipRadius(d); })
		.style("fill", '#BBB')
		.style("stroke",'#BBB')
		.on('mouseover', function(d) {
			virusTooltip.show(d, this);
		})
		.on('click', function(d) {
			if ((typeof d.db != "undefined") && (d.db == "GISAID") && (typeof d.accession != "undefined")) {
				var url = "http://gisaid.org/EPI/"+d.accession;
				console.log("opening url "+url);
				var win = window.open(url, '_blank');
  				win.focus();
  			}	
  		})		
		.on('mouseout', virusTooltip.hide);

	if (tips.length<100){
		var labels = treeplot.selectAll(".label")
			.data(tips)
			.enter()
			.append("text")
			.attr("class","label")
			.attr("x", function(d) { return d.x+10; })
			.attr("y", function(d) { return d.y+4; })
			.text(function(d) { return d.strain;});
		}

	var mutations = treeplot.selectAll(".muts")
		.data(nodes)
		.enter()
		.append("text")
		.attr("class", "")
		.attr("x", function(d) {
			return d.x - 6;
		})
		.attr("y", function(d) {
			return d.y - 3;
		})
		.style("text-anchor", "end")
		.text(function (d) {
			return d.aa_muts.replace(/,/g, ', ');
		});


	d3.select("#reset")
		.on("click", function(d) {
			var dMin = d3.min(xValues),
				dMax = d3.max(xValues),
				lMin = d3.min(yValues),
				lMax = d3.max(yValues);
			rescale(dMin, dMax, lMin, lMax);
		});

	function rescale(dMin, dMax, lMin, lMax) {
		var speed = 1500;
		xScale.domain([dMin,dMax]);
		yScale.domain([lMin,lMax]);

		nodes.forEach(function (d) {
			d.x = xScale(d.xvalue);
			d.y = yScale(d.yvalue);
		});

		treeplot.selectAll(".tip").data(tips)
			.transition().duration(speed)
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });

		treeplot.selectAll(".internal").data(internals)
			.transition().duration(speed)
			.attr("x", function(d) {
				if (typeof d.frequency != "undefined") {
					return d.x - 5*Math.sqrt(d.frequency) - 0.5;
				}
				else {
					return d.x - 1;
				}
			})
			.attr("y", function(d) {
				if (typeof d.frequency != "undefined") {
					return d.y - 5*Math.sqrt(d.frequency) - 0.5;
				}
				else {
					return d.y - 1;
				}
			});

		treeplot.selectAll(".link").data(links)
			.transition().duration(speed)
			.attr("points", function(d) {
				return (d.source.x).toString() + "," + d.source.y.toString() + " "
				+ (d.source.x).toString() + "," + d.target.y.toString() + " "
				+ (d.target.x).toString() + "," + d.target.y.toString()
			});
			
		treeplot.selectAll(".annotation").data(clades)
			.transition().duration(speed)
			.attr("x", function(d) {
				return xScale(d[1]) - 6;
			})
			.attr("y", function(d) {
				return yScale(d[2]) - 6;
			});			
	}	

	d3.select(window).on('resize', resize); 
	
	function resize() {
	
		var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
		var width = containerWidth,
			height = treePlotHeight(containerWidth);
			
		d3.select("#treeplot")
			.attr("width", width)
			.attr("height", height);			
			
		xScale.range([10, width-10]);
		yScale.range([10, height-10]);
		
		nodes.forEach(function (d) {
			d.x = xScale(d.xvalue);
			d.y = yScale(d.yvalue);
		});		
		
		treeplot.selectAll(".tip").data(tips)
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });

		treeplot.selectAll(".link").data(links)
			.attr("points", function(d) {
				return (d.source.x).toString() + "," + d.source.y.toString() + " "
				+ (d.source.x).toString() + "," + d.target.y.toString() + " "
				+ (d.target.x).toString() + "," + d.target.y.toString()
			});
			
		treeplot.selectAll(".annotation").data(clades)
			.attr("x", function(d) {
				return xScale(d[1]) - 6;
			})
			.attr("y", function(d) {
				return yScale(d[2]) - 6;
			});

	}
	
	function colorByGenotype() {
		var positions_string = document.getElementById("gt-color").value.split(',');
		var positions_list = []
		positions_string.map(function(d) {
			val = parseInt(d)-1;
			if (!isNaN(val)) {
				if (val < 551) {
					positions_list.push(val);
				}
			}
		});
		console.log(positions_list);
		if (positions_list.length > 0) {
			colorBy = "genotype";
		}
		else {
			positions_list = [1];
		}
		colorByGenotypePosition(positions_list);
	}

	function colorByGenotypePosition (positions) {
		var gts = nodes.map(function (d) {var tmp = [];
											for (var i=0; i<positions.length; i++){
												var aa = cladeToSeq[d.clade];
												tmp[tmp.length] = (positions[i]+1)+aa[positions[i]];
											}
											d.color_gt = tmp.join(" / "); 
											return d.color_gt;});
		var unique_gts = d3.set(gts).values();
		var gt_counts = {};
		for (var i=0; i<unique_gts.length; i++){gt_counts[unique_gts[i]]=0;}
		gts.forEach(function (d) {gt_counts[d]+=1;});
		var filtered_gts = unique_gts.filter(function (d) {return gt_counts[d]>=10;});
		filtered_gts.sort(function (a,b){var res;
			if (gt_counts[a]>gt_counts[b]){ res=-1;}
			else if (gt_counts[a]<gt_counts[b]){ res=1;}
			else {res=0;}
			return res;});
		console.log("genotypes passed filtering:"+filtered_gts);
		colorScale = d3.scale.ordinal()
			.domain(filtered_gts)
			.range(genotypeColors);			
 		treeplot.selectAll(".link")
			.style("stroke", function(d) {
				var col = colorScale(d.target.color_gt);
				return branchStrokeColor(col);
			});
 		treeplot.selectAll(".tip")
			.style("fill", function(d) {
				var col = colorScale(d.color_gt);
				return tipFillColor(col);
			})
			.style("stroke", function(d) {
				var col = colorScale(d.color_gt);
				return tipStrokeColor(col);
			});
		if (typeof tree_legend != undefined){
			removeLegend();
		}
		tree_legend = makeLegend();
	}
	
	function onSelect(tip) {
		d3.select("#"+(tip.strain).replace(/\//g, ""))
			.call(function(d) {
				virusTooltip.show(tip, d[0][0]);
			});
	}

	var mc = autocomplete(document.getElementById('search'))
		.keys(tips)
		.dataField("strain")
		.placeHolder("search strains...")
		.width(800)
		.height(500)
		.onSelected(onSelect)
		.render();



	d3.select("#gt-color")
		.on("keyup", colorByGenotype);

	colorByGenotype();
	tree_legend = makeLegend();

});

d3.json(file_prefix + "meta.json", function(error, json) {
	if (error) return console.warn(error);
	d3.select("#updated").text(json['updated']);
	commit_id = json['commit'];
	short_id = commit_id.substring(0, 6);	
	d3.select("#commit")
		.append("a")
		.attr("href", "http://github.com/blab/nextflu/commit/" + commit_id)
		.text(short_id);

});

d3.json(file_prefix + "sequences.json", function(error, json) {
	if (error) return console.warn(error);
	cladeToSeq=json;
});

