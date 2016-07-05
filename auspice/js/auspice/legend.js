var legendRectSize = 15;
var legendSpacing = 4;
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
		if ((colorBy=='lbi')||(colorBy=='date')||(colorBy=='dfreq')||(colorBy=='HI_dist')||(colorBy=='cHI')){
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

function removeLegend(){
	legend.selectAll('.legend')
  .remove();
}

