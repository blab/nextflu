function parseSequences(){
    var lines = document.getElementById('seqinput').value.split('\n');
    var seqs = {};
    var closest_nodes = {};
    var current_seq_name = "";
    var current_seq = "";
    for (var li=0; li<lines.length; li++){
        if (lines[li][0]=='>'){
            if (current_seq.length){
                current_seq_name += " seq "+Object.keys(seqs).length;
                seqs[current_seq_name]=current_seq;
             }
            current_seq_name = lines[li].substring(1,lines[li].length);
            current_seq = "";
        }else{
            current_seq += lines[li].toUpperCase().replace(/[^ACGTWRN]/g,"");
        }
    }
    if (current_seq.length){
        current_seq_name += " seq "+Object.keys(seqs).length;
        seqs[current_seq_name]=current_seq;
    }
    for (current_seq_name in seqs){
        var tmpclade = locateSequence(current_seq_name, seqs[current_seq_name]);
        if (typeof closest_nodes[tmpclade]=="undefined"){closest_nodes[tmpclade]=[current_seq_name];}
        else{closest_nodes[tmpclade].push(current_seq_name);}
    }
    for (var tmpclade in closest_nodes){
        markInTree(tmpclade);
    }
}

function locateSequence(name, seq){
    var mutations, olap_start, olap_end;
    console.log('Provided sequence: '+ name +': ' + seq.substring(0,20)+'....');
    tmp = alignToRoot(seq);
    olap_start=tmp[0]; olap_end=tmp[1]; mutations=tmp[2];
    console.log("start, end:", olap_start, olap_end);
    console.log("mutations:", mutations);
    var bestClade = findClosestClade(mutations);
    return bestClade;
}

function findClosestClade(mutations){
    var bestClade=-1, bestScore=0;
    var tmpScore=0;
    for (clade in cladeToSeq){
        tmpScore=0;
        for (mut in mutations){
            if (stateAtPosition(clade, 'nuc', mut)==mutations[mut]){
                tmpScore++;
            }
        }
        if (clade!="root") {
            tmpScore -= 0.5*Object.keys(cladeToSeq[clade]['nuc']).length;
        }
        if (tmpScore>bestScore){
            bestScore=tmpScore;
            bestClade=clade;
        }
    }
    console.log("best match:",bestClade);
    return bestClade;
}

function markInTree(clade){
    treeplot.selectAll('.tip').filter(function(d){return d.clade==clade;})
        .attr("r", function(d){console.log(d.strain); return tipRadius*2.7;})
        .style("visibility", 'visible')
        .style("fill", function (t) {
          return d3.rgb(tipFillColor(t)).brighter();
        });
}


function alignToRoot(seq){
    var rootSeq = cladeToSeq["root"]["nuc"];
    var shift = 0;
    var max_score = 0.0, max_shift;

    for(shift=0; shift<seq.length-30;shift++){
        var tmp_score = 0;
        var olaplen=Math.min(seq.length-shift, rootSeq.length);
        for (var pos=0; pos<olaplen; pos++){
            if (seq[pos+shift]==rootSeq[pos]){
                tmp_score++;
            }
        }
        tmp_score*=1.0/olaplen;
        if (tmp_score>max_score){
            max_score=tmp_score;
            max_shift=-shift;
        }
    }

    for(shift=0; shift<rootSeq.length-30;shift++){
        var tmp_score = 0;
        var olaplen=Math.min(rootSeq.length-shift, seq.length);
        for (var pos=0; pos<olaplen; pos++){
            if (seq[pos]==rootSeq[shift+pos]){
                tmp_score++;
            }
        }
        tmp_score*=1.0/olaplen;
        if (tmp_score>max_score){
            max_score=tmp_score;
            max_shift=shift;
        }
    }
    console.log("best shift: ",max_shift, " score: ",max_score);
    var mutations = {};
    if (max_shift<0){
        var olaplen=Math.min(seq.length-max_shift, rootSeq.length);
        var olap_start = 0;
        var olap_end = olaplen;
    }else{
        var olaplen=Math.min(rootSeq.length-max_shift, seq.length);
        var olap_start = max_shift;
        var olap_end = max_shift+olaplen;
    }
    for (var pos=olap_start; pos<olap_end; pos++){
        if (rootSeq[pos]!=seq[pos-max_shift]){
            mutations[pos]=seq[pos-max_shift];
        }
    }
    return [olap_start, olap_end, mutations];
}



d3.select('#seqinputsubmit').on('click', parseSequences); 
