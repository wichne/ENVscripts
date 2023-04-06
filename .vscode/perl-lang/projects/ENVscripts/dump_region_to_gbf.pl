{"version":5,"vars":[{"line":2,"containerName":"","kind":2,"name":"strict"},{"name":"%ENV","kind":13,"containerName":null,"line":3},{"line":3,"kind":12,"name":"ENVSCRIPTS"},{"kind":2,"name":"lib","line":3,"containerName":""},{"name":"ENV","kind":2,"containerName":"","line":4},{"line":5,"containerName":"","kind":2,"name":"DBI"},{"containerName":"Getopt","line":6,"name":"Std","kind":2},{"containerName":"Bio","line":7,"name":"SeqIO","kind":2},{"kind":2,"name":"Seq","line":8,"containerName":"Bio"},{"name":"LocationI","kind":2,"containerName":"Bio","line":9},{"containerName":null,"name":"%opt","localvar":"my","line":11,"kind":13,"definition":"my"},{"name":"getopts","kind":12,"line":12},{"kind":13,"name":"$opt","line":12,"containerName":null},{"definition":"my","kind":13,"line":13,"localvar":"my","name":"$db","containerName":null},{"line":13,"containerName":null,"kind":13,"name":"%opt"},{"containerName":null,"name":"$seqid","line":14,"localvar":"my","definition":"my","kind":13},{"kind":13,"name":"%opt","line":14,"containerName":null},{"localvar":"my","line":15,"kind":13,"definition":"my","containerName":null,"name":"$outfile"},{"name":"%opt","kind":13,"containerName":null,"line":15},{"line":16,"localvar":"my","definition":"my","kind":13,"containerName":null,"name":"$locoord"},{"containerName":null,"line":16,"name":"%opt","kind":13},{"containerName":null,"name":"$hicoord","localvar":"my","line":17,"definition":"my","kind":13},{"containerName":null,"line":17,"name":"%opt","kind":13},{"containerName":null,"name":"$host","line":19,"localvar":"my","kind":13,"definition":"my"},{"kind":13,"name":"%ENV","line":19,"containerName":null},{"line":19,"kind":12,"name":"DBSERVER"},{"kind":13,"name":"%ENV","line":19,"containerName":null},{"kind":12,"name":"DBSERVER","line":19},{"definition":"my","kind":13,"localvar":"my","line":20,"name":"$dbh","containerName":null},{"line":20,"name":"DBI","kind":12},{"kind":12,"name":"connect","line":20,"containerName":"main::"},{"line":22,"localvar":"my","definition":"my","kind":13,"containerName":null,"name":"$OUT"},{"line":23,"containerName":null,"kind":13,"name":"%outfile"},{"name":"$o","containerName":null,"definition":"my","kind":13,"line":24,"localvar":"my"},{"containerName":null,"line":24,"name":"$outfile","kind":13},{"containerName":null,"line":25,"name":"$OUT","kind":13},{"containerName":"SeqIO","line":25,"name":"Bio","kind":12},{"line":25,"containerName":"main::","kind":12,"name":"new"},{"kind":13,"name":"%outfile","line":28,"containerName":null},{"kind":13,"name":"$OUT","line":29,"containerName":null},{"kind":12,"name":"Bio","line":29,"containerName":"SeqIO"},{"name":"new","kind":12,"containerName":"main::","line":29},{"definition":"my","kind":13,"localvar":"my","line":33,"name":"$SeqObj","containerName":null},{"kind":12,"name":"seq_id_to_SeqObj","line":33},{"name":"$dbh","kind":13,"containerName":null,"line":33},{"name":"$seqid","kind":13,"containerName":null,"line":33},{"name":"add_features_to_SeqObj","kind":12,"line":34},{"containerName":null,"line":34,"name":"$dbh","kind":13},{"name":"$seqid","kind":13,"containerName":null,"line":34},{"containerName":null,"line":34,"name":"$SeqObj","kind":13},{"definition":"my","kind":13,"localvar":"my","line":36,"name":"$subseq","containerName":null},{"kind":13,"name":"$SeqObj","line":36,"containerName":null},{"containerName":"main::","line":36,"name":"subseq","kind":12},{"line":36,"containerName":null,"kind":13,"name":"$locoord"},{"line":36,"containerName":null,"kind":13,"name":"$hicoord"},{"containerName":null,"name":"$ssobj","localvar":"my","line":37,"kind":13,"definition":"my"},{"kind":12,"name":"Bio","line":37,"containerName":"Seq"},{"containerName":"main::","line":37,"name":"new","kind":12},{"kind":13,"name":"$SeqObj","line":37,"containerName":null},{"kind":12,"name":"display_id","line":37,"containerName":"main::"},{"kind":13,"name":"$subseq","line":37,"containerName":null},{"containerName":null,"name":"@features","line":38,"localvar":"my","definition":"my","kind":13},{"kind":13,"name":"$SeqObj","line":38,"containerName":null},{"line":38,"containerName":"main::","kind":12,"name":"get_SeqFeatures"},{"name":"$featObj","containerName":null,"kind":13,"definition":"my","localvar":"my","line":39},{"name":"@features","kind":13,"containerName":null,"line":39},{"containerName":null,"name":"$start","localvar":"my","line":40,"kind":13,"definition":"my"},{"containerName":null,"line":40,"name":"$featObj","kind":13},{"containerName":"main::","line":40,"name":"start","kind":12},{"containerName":null,"name":"$end","line":41,"localvar":"my","definition":"my","kind":13},{"containerName":null,"line":41,"name":"$featObj","kind":13},{"kind":12,"name":"end","line":41,"containerName":"main::"},{"name":"$locoord","kind":13,"containerName":null,"line":42},{"kind":13,"name":"$end","line":42,"containerName":null},{"name":"$hicoord","kind":13,"containerName":null,"line":42},{"name":"$end","kind":13,"containerName":null,"line":42},{"containerName":null,"line":43,"name":"$locoord","kind":13},{"containerName":null,"line":43,"name":"$start","kind":13},{"containerName":null,"line":43,"name":"$hicoord","kind":13},{"kind":13,"name":"$start","line":43,"containerName":null},{"name":"$start","kind":13,"containerName":null,"line":44},{"kind":13,"name":"$locoord","line":44,"containerName":null},{"kind":13,"name":"$end","line":44,"containerName":null},{"containerName":null,"line":44,"name":"%hicoord","kind":13},{"name":"$thisFeat","containerName":null,"definition":"my","kind":13,"localvar":"my","line":45},{"line":45,"containerName":null,"kind":13,"name":"$featObj"},{"kind":13,"name":"$thisFeat","line":46,"containerName":null},{"kind":12,"name":"start","line":46,"containerName":"main::"},{"kind":13,"name":"$featObj","line":46,"containerName":null},{"name":"start","kind":12,"containerName":"main::","line":46},{"line":46,"containerName":null,"kind":13,"name":"$locoord"},{"containerName":null,"line":47,"name":"$thisFeat","kind":13},{"name":"start","kind":12,"containerName":"main::","line":47},{"containerName":null,"line":48,"name":"$thisFeat","kind":13},{"kind":12,"name":"start","line":48,"containerName":"main::"},{"line":50,"containerName":null,"kind":13,"name":"$featObj"},{"containerName":"main::","line":50,"name":"end","kind":12},{"kind":13,"name":"%hicoord","line":50,"containerName":null},{"line":51,"containerName":null,"kind":13,"name":"$thisFeat"},{"containerName":"main::","line":51,"name":"end","kind":12},{"name":"$hicoord","kind":13,"containerName":null,"line":51},{"containerName":null,"line":51,"name":"%locoord","kind":13},{"name":"$thisFeat","kind":13,"containerName":null,"line":53},{"containerName":"main::","line":53,"name":"end","kind":12},{"name":"$featObj","kind":13,"containerName":null,"line":53},{"name":"end","kind":12,"containerName":"main::","line":53},{"line":53,"containerName":null,"kind":13,"name":"$locoord"},{"containerName":null,"line":55,"name":"$thisFeat","kind":13},{"containerName":"main::","line":55,"name":"primary_tag","kind":12},{"line":56,"containerName":null,"kind":13,"name":"$thisFeat"},{"name":"start","kind":12,"containerName":"main::","line":56},{"containerName":null,"line":57,"name":"$thisFeat","kind":13},{"line":57,"containerName":"main::","kind":12,"name":"end"},{"containerName":null,"line":57,"name":"$hicoord","kind":13},{"kind":13,"name":"$locoord","line":57,"containerName":null},{"kind":13,"name":"$ssobj","line":59,"containerName":null},{"containerName":"main::","line":59,"name":"add_SeqFeature","kind":12},{"line":59,"containerName":null,"kind":13,"name":"$thisFeat"},{"containerName":null,"line":63,"name":"$OUT","kind":13},{"name":"write_seq","kind":12,"containerName":"main::","line":63},{"name":"$ssobj","kind":13,"containerName":null,"line":63}]}