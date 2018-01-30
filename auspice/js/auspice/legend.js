var legendRectSize = 15;
var legendSpacing = 4;
var map_features;
function makeLegend(){

	d3.select("#legend-title").html(function(d){
		if (colorBy == "ep") {
			return "Epitope mutations";
		}
		if (colorBy == "ne") {
			return "Non-epitope mutations";
		}
		if (colorBy == "rb") {
			return "Receptor binding mutations";
		}
		if (colorBy == "lbi") {
			return "Local branching index";
		}
		if (colorBy == "region") {
			return "Region";
		}
		if (colorBy == "genotype") {
			return "Genotype";
		}
		if (colorBy == "date") {
			return "Date";
		}
        if (colorBy == "cHI") {
            return "log<sub>2</sub> titer distance from root";
        }
        if (colorBy == "HI_dist") {
            return "log<sub>2</sub> titer distance from "+focusNode.strain;
        }
		if (colorBy == "dfreq") {
			var tmp_nmonth = Math.round(12*dfreq_dn*time_step);
			var tmp_text = "Freq. change ("+tmp_nmonth+" month";
			if (tmp_nmonth>1){
				tmp_text+='s';
			}
			return tmp_text+')';
		}
		if (colorBy == "fitness") {
			return "Relative fitness";
		}
	});

	if (colorBy == "division"){
		make_map();
	} else {
		make_panels();
	}

}

function removeLegend(){
	legend.selectAll('.legend').remove();
	legend.selectAll('.map_feature').remove();
}

function make_panels(){

	// construct a dictionary that maps a legend entry to the preceding interval
	var lower_bound = {}, upper_bound = {};
	lower_bound[colorScale.domain()[0]] = -100000000;
    upper_bound[colorScale.domain()[0]] = colorScale.domain()[0];
	for (var i=1; i<colorScale.domain().length; i++){
		lower_bound[colorScale.domain()[i]]=colorScale.domain()[i-1];
        upper_bound[colorScale.domain()[i]]=colorScale.domain()[i];
	}
    upper_bound[colorScale.domain()[colorScale.domain().length-1]]=10000000;
	// function that equates a tip and a legend element
	// exact match is required for categorical qunantities such as genotypes, regions
	// continuous variables need to fall into the interal (lower_bound[leg], leg]
	var legend_match = function(leg, tip){
		if ((colorBy=='glyc')||(colorBy=='age')||(colorBy=='lbi')||(colorBy=='date')||(colorBy=='dfreq')||(colorBy=='HI_dist')||(colorBy=='cHI')){
			return (tip.coloring<=upper_bound[leg])&&(tip.coloring>lower_bound[leg]);
		}else{
			return tip.coloring==leg;
		}
	}

	var count = colorScale.domain().length;
	var stack = Math.ceil(count / 2);
	d3.select("#legend")
		.attr("height", stack * (legendRectSize + legendSpacing) + legendSpacing);

	var tmp_leg = legend.selectAll(".legend")
	.data(colorScale.domain())
	.enter().append('g')
	.attr('class', 'legend')
	.attr('transform', function(d, i) {
		var fromRight = Math.floor(i / stack);
		var fromTop = i % stack;
		var horz = fromRight * 145 + 5;
		var vert = fromTop * (legendRectSize + legendSpacing) + 5;
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
   		return d3.rgb(col).toString();
 	})
   .on('mouseover', function(leg){
    	treeplot.selectAll(".tip") //highlight all tips corresponding to legend
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius(d)*1.7;})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t)).brighter();
            });
		})
  	.on('mouseout', function(leg){
    	treeplot.selectAll(".tip") //undo highlight
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius(d);})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t));
            });
	    });

    tmp_leg.append('text')
    .attr('x', legendRectSize + legendSpacing + 5)
    .attr('y', legendRectSize - legendSpacing)
    .text(function(d) {
		var label = d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ').replace(/([a-z]+)_([a-z]+)/g, function(_, a, b) { return a.toTitleCase().concat(' ', b.toTitleCase()); }).replace(/^([a-z]+)$/, function(_, a) { return a.toTitleCase(); });
		label = label.replace(/([A-Za-z]+)_([A-Za-z]+)/g, function(_, a, b) { return a.toTitleCase().concat(' ', b.toTitleCase()); }).replace(/^([a-z]+)$/, function(_, a) { return a.toTitleCase(); });
		label = label.replace(/^Usa/, 'USA').replace(/^Usvi/, 'USVI');
    if (colorBy == "dfreq") {
        label += "\u00D7";
    }
    return label;
    })
   .on('mouseover', function(leg){
    	treeplot.selectAll(".tip")
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius(d)*1.7;})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t)).brighter();
            });
		})
  	.on('mouseout', function(leg){
    	treeplot.selectAll(".tip")
            .filter(function (d){return legend_match(leg, d);})
            .attr("r", function(d){return tipRadius(d);})
            .style("fill", function (t) {
              return d3.rgb(tipFillColor(t));
            });
	    });
	return tmp_leg;

}

function patch_division_name(d){
  var tmp = d.properties.NAME_2;
  if (tmp == null){
    tmp= d.properties.NAME_1;
  }
  if (tmp=='?'||tmp==null){
    tmp = d.properties.ISO;
    console.log('Falling back on ISO: '+ d.properties.NAME_2+ ' ' + d.properties.NAME_1 + ' ' + d.properties.ISO + ' ' +d.id);
  }
  return tmp.replace(' ','');
}

function patch_color(d) {
	return divisionColorScale(patch_division_name(d));
}

function make_map(){

    var width = 380,
        height = 320,
        active = d3.select(null);

    console.log('enter map');
    var projection = d3.geo.robinson()
            .scale(3000)
            .translate([width, height * 1.5])

    var path = d3.geo.path()
        .projection(projection);

    var svg = d3.select("#legend");
    svg.call(mapTooltip);

    svg
        .attr("width", width)
        .attr("height", height);

    svg.append("rect")
        .attr("class", "map_background")
        .attr("width", width)
        .attr("height", height);

    var g = svg.append("g")
        .style("stroke-width", "2px");

    d3.json("/data/ebola_map.json", function(error, locations) {
        var locationData = topojson.feature(locations, locations.objects.ebov).features;

				console.log(locationData);

        map_features = g.selectAll(".map_feature")
            .data(locationData)
            .enter().append("path")
            .style("fill", patch_color)
            .attr("d", path)
            .attr("class", "map_feature")
            .on("mouseover",mouseOverMap)
            .on("mouseout",mouseOutMap);

        g.attr("transform", "translate(260,180) scale(" + 0.65 + ")");

      g.append("path")
          .datum(topojson.mesh(locations, locations.objects.ebov, function(a, b) { return a !== b; }))
          .attr("class", "map_mesh")
          .attr("d", path);

    });
}

function match_division(map_division, tip){
  var tmp = patch_division_name(map_division);
  return tmp==tip.attr.division;
}

function mouseOverMap(division){
  mapTooltip.show(division);
  treeplot.selectAll(".tip")
		.filter(function (d){ return match_division(division, d);})
    .attr("r", function(d){return tipRadius(d)*1.7;})
		.style("fill", function (t) {
			return d3.rgb(tipFillColor(t)).brighter();
		});
	legend.selectAll('.map_feature')
		.filter(function (m) { return patch_division_name(m) == patch_division_name(division);})
		.style("fill", function(m) {
			return d3.rgb(patch_color(division)).brighter();
		});
}

function mouseOutMap(division){
    mapTooltip.hide(division);
    treeplot.selectAll(".tip")
            .filter(function (d){ return match_division(division, d);})
                .attr("r", function(d){return tipRadius(d);})
								.style("fill", function (t) {
									return d3.rgb(tipFillColor(t));
								});
	legend.selectAll('.map_feature')
		.filter(function (m) { return patch_division_name(m) == patch_division_name(division);})
		.style("fill", function(m) {
			return d3.rgb(patch_color(division));
		});
}
