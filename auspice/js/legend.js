var legendRectSize = 15;
var legendSpacing = 4;
function makeLegend(){
	
	d3.select("#legend-title").text(function(d){
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
		if (colorBy == "country") {
			return "Country";
		}
        if (colorBy == "host") {
            return "Host";
        }
		if (colorBy == "genotype") {
			return "Genotype";
		}
		if (colorBy == "date") {
			return "Date";
		}		
		if (colorBy == "dfreq") {
			var tmp_nmonth = Math.round(12*dfreq_dn*time_step);
			var tmp_text = "Freq. change ("+tmp_nmonth+" month";
			if (tmp_nmonth>1){
				tmp_text+='s';
			}
			return tmp_text+')';
		}
});
 console.log(colorScale.domain());
 var dd;
 if (colorBy=='date'){
  dd =  colorScale.domain()[1]-colorScale.domain()[0];
 }else{
  dd=0;
 }
 var stack = Math.ceil(colorScale.domain().filter(function(d){return typeof d != "undefined";}).length/2);
  var tmp_leg = legend.selectAll(".legend")
  .data(colorScale.domain().filter(function(d){return typeof d != "undefined";}))
  .enter().append('g')
  .attr('class', 'legend')
  .attr('transform', function(d, i) {
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
   return d3.rgb(col).toString();
 })
  .on('mouseover', function(leg){
    treeplot.selectAll(".tip")
            .filter(function (d){
              if (colorBy=='date') {return d.coloring<leg && d.coloring>=leg-dd;}
              else if (colorBy=='host' || colorBy=="country") {return d.coloring==leg;}
            })
                .attr("r", function(d){return tipRadius*1.7;})
                .style("fill", function (t) {
                  return d3.rgb(tipFillColor(t)).brighter();
                });
  }) 
  .on('mouseout', function(leg){
    treeplot.selectAll(".tip")
            .filter(function (d){
              if (colorBy=='date') {return d.coloring<leg && d.coloring>=leg-dd;}
              else if (colorBy=='host' || colorBy=="country") {return d.coloring==leg;}
            })
                .attr("r", function(d){return tipRadius;})
                .style("fill", function (t) {
                  return d3.rgb(tipFillColor(t));
                });
    });


  tmp_leg.append('text')
  .attr('x', legendRectSize + legendSpacing + 5)
  .attr('y', legendRectSize - legendSpacing)
  .text(function(d) {
    if (colorBy=="country"){
        console.log(d);
       return codeToCountry[d];
    }else{
       return d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');
    }
    })
  .on('mouseover', function(leg){
    treeplot.selectAll(".tip")
            .filter(function (d){
              if (colorBy=='date') {return d.coloring<leg && d.coloring>=leg-dd;}
              else if (colorBy=='host' || colorBy=="country") {return d.coloring==leg;}
            })
                .attr("r", function(d){return tipRadius*1.7;})
                .style("fill", function (t) {
                  return d3.rgb(tipFillColor(t)).brighter();
                });
  }) 
  .on('mouseout', function(leg){
    treeplot.selectAll(".tip")
            .filter(function (d){
              if (colorBy=='date') {return d.coloring<leg && d.coloring>=leg-dd;}
              else if (colorBy=='host' || colorBy=="country") {return d.coloring==leg;}
            })
                .attr("r", function(d){return tipRadius;})
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

