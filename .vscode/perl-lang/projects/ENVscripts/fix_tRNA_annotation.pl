{"version":5,"vars":[{"containerName":null,"line":2,"kind":13,"name":"%ENV"},{"kind":12,"line":2,"name":"SCRIPTS"},{"kind":2,"line":2,"containerName":"","name":"lib"},{"containerName":"","line":3,"kind":2,"name":"ENV"},{"containerName":"","line":4,"kind":2,"name":"strict"},{"kind":2,"line":5,"containerName":"Getopt","name":"Std"},{"kind":2,"line":6,"containerName":"Bio::SeqFeature","name":"Generic"},{"name":"DBI","line":7,"kind":2,"containerName":""},{"name":"%args","localvar":"my","line":9,"definition":"my","kind":13,"containerName":null},{"line":10,"kind":12,"name":"getopts"},{"name":"$args","containerName":null,"kind":13,"line":10},{"containerName":null,"kind":13,"name":"$DEBUG","definition":"my","localvar":"my","line":11},{"name":"$input_file","definition":"my","localvar":"my","line":12,"containerName":null,"kind":13},{"line":12,"kind":13,"containerName":null,"name":"%args"},{"kind":13,"containerName":null,"name":"$db","line":13,"localvar":"my","definition":"my"},{"containerName":null,"kind":13,"line":13,"name":"%args"},{"kind":13,"containerName":null,"line":14,"localvar":"my","definition":"my","name":"$pwd"},{"line":14,"kind":13,"containerName":null,"name":"%args"},{"name":"$user","line":15,"localvar":"my","definition":"my","kind":13,"containerName":null},{"kind":13,"line":15,"containerName":null,"name":"%args"},{"line":15,"kind":13,"containerName":null,"name":"%args"},{"name":"%ENV","containerName":null,"line":15,"kind":13},{"line":15,"kind":12,"name":"USER"},{"containerName":null,"kind":13,"definition":"my","localvar":"my","line":16,"name":"$host"},{"name":"%ENV","line":16,"kind":13,"containerName":null},{"name":"DBSERVER","kind":12,"line":16},{"containerName":null,"kind":13,"line":16,"name":"%ENV"},{"name":"DBSERVER","kind":12,"line":16},{"definition":"my","localvar":"my","line":18,"name":"$dbh","containerName":null,"kind":13},{"line":18,"kind":12,"name":"DBI"},{"name":"connect","line":18,"kind":12,"containerName":"main::"},{"name":"$user","kind":13,"line":18,"containerName":null},{"line":18,"kind":13,"containerName":null,"name":"$pwd"},{"localvar":"my","line":20,"definition":"my","name":"$IN","kind":13,"containerName":null},{"line":20,"kind":13,"containerName":null,"name":"$input_file"},{"line":22,"localvar":"my","definition":"my","name":"%AA","kind":13,"containerName":null},{"kind":13,"containerName":null,"line":23,"localvar":"my","definition":"my","name":"$line"},{"name":"%IN","line":23,"kind":13,"containerName":null},{"kind":13,"line":24,"containerName":null,"name":"$line"},{"definition":"my","localvar":"my","line":25,"name":"$seq_acc","containerName":null,"kind":13},{"name":"$index","kind":13,"line":26,"containerName":null},{"name":"$end5","containerName":null,"line":27,"kind":13},{"containerName":null,"kind":13,"line":28,"name":"$end3"},{"name":"$aa","line":29,"kind":13,"containerName":null},{"line":30,"kind":13,"containerName":null,"name":"$acodon"},{"name":"$intron_end5","containerName":null,"kind":13,"line":31},{"name":"$intron_end3","kind":13,"line":32,"containerName":null},{"containerName":null,"line":33,"kind":13,"name":"$score"},{"name":"$line","line":33,"kind":13,"containerName":null},{"kind":13,"line":34,"containerName":null,"name":"%index"},{"kind":12,"line":34,"name":"next"},{"containerName":null,"kind":13,"name":"@accs","definition":"my","line":35,"localvar":"my"},{"containerName":null,"kind":13,"line":35,"name":"$seq_acc"},{"containerName":null,"kind":13,"name":"$sid","definition":"my","localvar":"my","line":36},{"containerName":null,"kind":13,"name":"$sacc","definition":"my","line":37,"localvar":"my"},{"name":"$ax","definition":"my","localvar":"my","line":38,"containerName":null,"kind":13},{"line":38,"kind":13,"containerName":null,"name":"@accs"},{"name":"$sid","containerName":null,"kind":13,"line":39},{"name":"get_seq_id_by_seq_accession","kind":12,"line":39},{"kind":13,"line":39,"containerName":null,"name":"$dbh"},{"name":"$ax","line":39,"kind":13,"containerName":null},{"name":"$sacc","containerName":null,"kind":13,"line":40},{"kind":13,"line":40,"containerName":null,"name":"$ax"},{"name":"%sid","containerName":null,"line":41,"kind":13},{"name":"last","kind":12,"line":41},{"kind":13,"line":43,"containerName":null,"name":"%sid"},{"line":45,"localvar":"my","definition":"my","name":"$start","kind":13,"containerName":null},{"kind":13,"line":45,"containerName":null,"name":"$end"},{"name":"$strand","kind":13,"line":45,"containerName":null},{"containerName":null,"line":45,"kind":13,"name":"$end5"},{"line":45,"kind":13,"containerName":null,"name":"$end3"},{"name":"$end5","containerName":null,"line":45,"kind":13},{"containerName":null,"line":45,"kind":13,"name":"$end3"},{"name":"$end3","containerName":null,"line":45,"kind":13},{"name":"$end5","containerName":null,"line":45,"kind":13},{"name":"%AA","containerName":null,"kind":13,"line":46},{"line":46,"kind":13,"containerName":null,"name":"$aa"},{"kind":13,"containerName":null,"localvar":"my","line":47,"definition":"my","name":"$id"},{"containerName":null,"line":47,"kind":13,"name":"$sacc"},{"containerName":null,"kind":13,"name":"$name","definition":"my","localvar":"my","line":48},{"containerName":null,"kind":13,"definition":"my","line":50,"localvar":"my","name":"$feature_q"},{"definition":"my","localvar":"my","line":51,"name":"$fid","containerName":null,"kind":13},{"name":"$dbh","containerName":null,"line":51,"kind":13},{"containerName":"main::","kind":12,"line":51,"name":"selectrow_array"},{"kind":13,"line":51,"containerName":null,"name":"$feature_q"},{"line":52,"kind":13,"containerName":null,"name":"%fid"},{"name":"$DEBUG","line":53,"kind":13,"containerName":null},{"name":"$annObj","localvar":"my","line":55,"definition":"my","kind":13,"containerName":null},{"name":"Bio","kind":12,"line":55,"containerName":"Annotation::Collection"},{"containerName":"main::","line":55,"kind":12,"name":"new"},{"name":"$annObj","kind":13,"line":56,"containerName":null},{"name":"add_Annotation","containerName":"main::","line":56,"kind":12},{"name":"Bio","containerName":"Annotation::SimpleValue","line":57,"kind":12},{"name":"new","containerName":"main::","line":57,"kind":12},{"name":"$name","containerName":null,"line":57,"kind":13},{"containerName":null,"kind":13,"line":58,"name":"$annObj"},{"name":"add_Annotation","containerName":"main::","line":58,"kind":12},{"name":"Bio","containerName":"Annotation::SimpleValue","line":59,"kind":12},{"name":"new","line":59,"kind":12,"containerName":"main::"},{"line":59,"kind":13,"containerName":null,"name":"$acodon"},{"name":"delete_feature_annotations","line":61,"kind":12},{"name":"$dbh","kind":13,"line":61,"containerName":null},{"containerName":null,"kind":13,"line":61,"name":"$fid"},{"kind":12,"line":62,"name":"load_feature_annotations"},{"kind":13,"line":62,"containerName":null,"name":"$dbh"},{"name":"$fid","containerName":null,"line":62,"kind":13},{"line":62,"kind":13,"containerName":null,"name":"$annObj"}]}