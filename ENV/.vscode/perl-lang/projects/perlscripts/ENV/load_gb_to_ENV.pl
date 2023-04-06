{"version":5,"vars":[{"kind":2,"containerName":"","line":2,"name":"DBI"},{"line":3,"containerName":"Bio","kind":2,"name":"SeqIO"},{"name":"Seq","line":4,"containerName":"Bio","kind":2},{"name":"Std","kind":2,"line":5,"containerName":"Getopt"},{"kind":2,"line":6,"containerName":"Data","name":"Dumper"},{"name":"strict","line":7,"containerName":"","kind":2},{"containerName":null,"line":8,"kind":13,"name":"%ENV"},{"name":"ENVSCRIPTS","kind":12,"line":8},{"name":"lib","containerName":"","line":8,"kind":2},{"name":"ENV","kind":2,"line":9,"containerName":""},{"name":"vars","line":10,"containerName":"","kind":2},{"definition":"my","name":"%args","localvar":"my","kind":13,"containerName":null,"line":12},{"line":13,"kind":12,"name":"getopts"},{"name":"$args","kind":13,"containerName":null,"line":13},{"name":"$filename","definition":"my","line":15,"containerName":null,"kind":13,"localvar":"my"},{"name":"%args","kind":13,"line":15,"containerName":null},{"name":"$format","definition":"my","line":16,"containerName":null,"localvar":"my","kind":13},{"containerName":null,"line":16,"kind":13,"name":"%args"},{"kind":13,"containerName":null,"line":16,"name":"%args"},{"line":17,"containerName":null,"localvar":"my","kind":13,"definition":"my","name":"$db"},{"name":"%args","line":17,"containerName":null,"kind":13},{"line":18,"containerName":null,"localvar":"my","kind":13,"name":"$pswd","definition":"my"},{"name":"%args","containerName":null,"line":18,"kind":13},{"kind":13,"localvar":"my","line":19,"containerName":null,"definition":"my","name":"$user"},{"name":"%args","containerName":null,"line":19,"kind":13},{"line":19,"containerName":null,"kind":13,"name":"%args"},{"name":"%ENV","kind":13,"containerName":null,"line":19},{"name":"$locus_prefix","definition":"my","kind":13,"localvar":"my","line":20,"containerName":null},{"line":20,"containerName":null,"kind":13,"name":"%args"},{"name":"$set_id","definition":"my","line":21,"containerName":null,"localvar":"my","kind":13},{"kind":13,"containerName":null,"line":21,"name":"%args"},{"kind":13,"localvar":"my","containerName":null,"line":22,"name":"$source","definition":"my"},{"line":22,"containerName":null,"kind":13,"name":"%args"},{"definition":"my","name":"$in","localvar":"my","kind":13,"containerName":null,"line":24},{"kind":12,"containerName":"SeqIO","line":24,"name":"Bio"},{"name":"new","kind":12,"line":24,"containerName":"main::"},{"name":"$filename","line":24,"containerName":null,"kind":13},{"line":25,"containerName":null,"kind":13,"name":"$format"},{"localvar":"my","kind":13,"containerName":null,"line":26,"name":"$host","definition":"my"},{"line":26,"containerName":null,"kind":13,"name":"%ENV"},{"line":26,"kind":12,"name":"DBSERVER"},{"kind":13,"containerName":null,"line":26,"name":"%ENV"},{"line":26,"kind":12,"name":"DBSERVER"},{"containerName":null,"line":27,"kind":13,"name":"$dbh"},{"name":"DBI","kind":12,"line":27},{"containerName":"main::","line":27,"kind":12,"name":"connect"},{"containerName":null,"line":28,"kind":13,"name":"$user"},{"line":28,"containerName":null,"kind":13,"name":"%pswd"},{"kind":13,"localvar":"my","line":30,"containerName":null,"definition":"my","name":"%DONE"},{"localvar":"my","kind":13,"line":32,"containerName":null,"definition":"my","name":"$seqrec"},{"containerName":null,"line":32,"kind":13,"name":"$in"},{"line":32,"containerName":"main::","kind":12,"name":"next_seq"},{"name":"@features","definition":"my","containerName":null,"line":33,"localvar":"my","kind":13},{"containerName":null,"line":33,"kind":13,"name":"$seqrec"},{"name":"get_SeqFeatures","kind":12,"line":33,"containerName":"main::"},{"definition":"my","name":"$dbseqobj","containerName":null,"line":37,"kind":13,"localvar":"my"},{"name":"$seq_id","kind":13,"containerName":null,"line":37},{"line":38,"containerName":null,"localvar":"my","kind":13,"definition":"my","name":"$seq_q"},{"name":"$seqrec","containerName":null,"line":40,"kind":13},{"name":"display_name","line":41,"containerName":"main::","kind":12},{"definition":"my","name":"$seq_q_r","line":44,"containerName":null,"localvar":"my","kind":13},{"kind":13,"line":44,"containerName":null,"name":"$dbh"},{"containerName":"main::","line":44,"kind":12,"name":"selectall_hashref"},{"kind":13,"line":44,"containerName":null,"name":"$seq_q"},{"line":45,"containerName":null,"kind":13,"name":"%seq_q_r"},{"name":"$seqrec","containerName":null,"line":45,"kind":13},{"kind":12,"line":45,"containerName":"main::","name":"display_name"},{"definition":"my","name":"$sid","containerName":null,"line":46,"localvar":"my","kind":13},{"name":"%seq_q_r","kind":13,"line":46,"containerName":null},{"kind":13,"containerName":null,"line":47,"name":"$seq_id"},{"line":47,"containerName":null,"kind":13,"name":"$sid"},{"name":"$dbseqobj","line":48,"containerName":null,"kind":13},{"line":48,"containerName":"Seq","kind":12,"name":"Bio"},{"name":"new","line":48,"containerName":"main::","kind":12},{"kind":13,"line":48,"containerName":null,"name":"%seq_q_r"},{"name":"%sid","kind":13,"containerName":null,"line":48},{"name":"$sid","line":49,"containerName":null,"kind":13},{"kind":13,"line":50,"containerName":null,"name":"%seq_q_r"},{"name":"%sid","kind":13,"line":50,"containerName":null},{"line":54,"containerName":null,"kind":13,"name":"%dbseqobj"},{"kind":13,"line":55,"containerName":null,"name":"$seqrec"},{"kind":12,"containerName":"main::","line":55,"name":"display_name"},{"name":"$seq_id","containerName":null,"line":57,"kind":13},{"line":57,"kind":12,"name":"load_sequence_SeqFeature"},{"line":57,"containerName":null,"kind":13,"name":"$dbh"},{"name":"$seqrec","kind":13,"line":57,"containerName":null},{"name":"$seqrec","kind":13,"containerName":null,"line":58},{"containerName":"main::","line":58,"kind":12,"name":"display_name"},{"line":61,"containerName":null,"kind":13,"name":"%set_id"},{"definition":"my","name":"$name","containerName":null,"line":63,"localvar":"my","kind":13},{"kind":13,"line":63,"containerName":null,"name":"$filename"},{"name":"$name","kind":13,"containerName":null,"line":64},{"line":65,"containerName":null,"kind":13,"name":"$set_id"},{"kind":12,"line":65,"name":"load_sequence_sets"},{"name":"$dbh","kind":13,"line":65,"containerName":null},{"name":"$name","containerName":null,"line":65,"kind":13},{"kind":13,"containerName":null,"line":65,"name":"$seqrec"},{"kind":12,"containerName":"main::","line":65,"name":"desc"},{"name":"%set_id","containerName":null,"line":66,"kind":13},{"name":"link_seq_to_set","kind":12,"line":69},{"kind":13,"containerName":null,"line":69,"name":"$dbh"},{"name":"$set_id","line":69,"containerName":null,"kind":13},{"containerName":null,"line":69,"kind":13,"name":"$seq_id"},{"name":"$dbseqobj","kind":13,"containerName":null,"line":71},{"name":"get_sequence_by_seq_id","kind":12,"line":71},{"name":"$dbh","kind":13,"line":71,"containerName":null},{"kind":13,"containerName":null,"line":71,"name":"$seq_id"},{"containerName":null,"line":81,"localvar":"my","kind":13,"definition":"my","name":"$feature"},{"name":"@features","line":81,"containerName":null,"kind":13},{"definition":"my","name":"$feat_type","kind":13,"localvar":"my","line":82,"containerName":null},{"name":"$feature","line":82,"containerName":null,"kind":13},{"name":"primary_tag","line":82,"containerName":"main::","kind":12},{"containerName":null,"line":84,"kind":13,"name":"%feat_type"},{"name":"next","line":84,"kind":12},{"name":"$feat_type","line":86,"containerName":null,"kind":13},{"name":"$feature","containerName":null,"line":86,"kind":13},{"line":86,"containerName":"main::","kind":12,"name":"has_tag"},{"kind":12,"line":86,"name":"next"},{"name":"%data","definition":"my","localvar":"my","kind":13,"line":89,"containerName":null},{"definition":"my","name":"$min","localvar":"my","kind":13,"containerName":null,"line":90},{"line":90,"containerName":null,"kind":13,"name":"$max"},{"line":90,"containerName":null,"kind":13,"name":"$strand"},{"name":"$feature","kind":13,"containerName":null,"line":90},{"name":"start","containerName":"main::","line":90,"kind":12},{"name":"$feature","kind":13,"containerName":null,"line":90},{"name":"end","kind":12,"line":90,"containerName":"main::"},{"line":90,"containerName":null,"kind":13,"name":"$feature"},{"name":"strand","line":90,"containerName":"main::","kind":12},{"containerName":null,"line":91,"kind":13,"name":"$strand"},{"name":"%strand","kind":13,"containerName":null,"line":91},{"kind":13,"containerName":null,"line":92,"name":"$min"},{"name":"$max","containerName":null,"line":92,"kind":13},{"name":"$strand","line":92,"containerName":null,"kind":13},{"line":92,"containerName":null,"kind":13,"name":"$min"},{"name":"$max","kind":13,"line":92,"containerName":null},{"kind":13,"containerName":null,"line":92,"name":"$min"},{"kind":13,"containerName":null,"line":92,"name":"$max"},{"kind":13,"containerName":null,"line":92,"name":"$max"},{"kind":13,"line":92,"containerName":null,"name":"$min"},{"line":94,"containerName":null,"localvar":"my","kind":13,"name":"$FT_coords","definition":"my"},{"name":"$feature","containerName":null,"line":94,"kind":13},{"kind":12,"containerName":"main::","line":94,"name":"location"},{"name":"to_FTstring","line":94,"containerName":"main::","kind":12},{"definition":"my","name":"@tags","localvar":"my","kind":13,"containerName":null,"line":97},{"name":"$feature","kind":13,"line":97,"containerName":null},{"kind":12,"line":97,"containerName":"main::","name":"get_all_tags"},{"kind":13,"localvar":"my","containerName":null,"line":98,"definition":"my","name":"$locus_tag"},{"containerName":null,"line":98,"kind":13,"name":"$max_partial"},{"kind":13,"containerName":null,"line":98,"name":"$min_partial"},{"name":"$tag","definition":"my","line":99,"containerName":null,"localvar":"my","kind":13},{"name":"@tags","line":99,"containerName":null,"kind":13},{"line":100,"containerName":null,"kind":13,"localvar":"my","name":"@values","definition":"my"},{"containerName":null,"line":100,"kind":13,"name":"$feature"},{"name":"get_tag_values","kind":12,"line":100,"containerName":"main::"},{"name":"$tag","kind":13,"line":100,"containerName":null},{"name":"%tag","kind":13,"containerName":null,"line":101},{"name":"$locus_tag","kind":13,"containerName":null,"line":101},{"containerName":null,"line":101,"kind":13,"name":"@values"},{"containerName":null,"line":102,"kind":13,"name":"%tag"},{"line":103,"containerName":null,"kind":13,"name":"$min_partial"},{"name":"@values","line":103,"containerName":null,"kind":13},{"kind":13,"line":104,"containerName":null,"name":"$max_partial"},{"name":"@values","kind":13,"line":104,"containerName":null},{"kind":13,"localvar":"my","containerName":null,"line":108,"definition":"my","name":"$feat_r"},{"kind":13,"line":109,"containerName":null,"name":"%locus_tag"},{"kind":13,"line":109,"containerName":null,"name":"$locus_tag"},{"kind":13,"line":109,"containerName":null,"name":"$feature"},{"containerName":"main::","line":109,"kind":12,"name":"primary_tag"},{"containerName":null,"line":113,"kind":13,"localvar":"my","definition":"my","name":"$acc_q"},{"kind":13,"containerName":null,"line":114,"name":"$feat_r"},{"name":"$dbh","containerName":null,"line":114,"kind":13},{"name":"selectcol_arrayref","containerName":"main::","line":114,"kind":12},{"name":"$acc_q","line":114,"containerName":null,"kind":13},{"name":"$feat_r","line":118,"containerName":null,"kind":13},{"line":118,"containerName":null,"kind":13,"name":"%feat_r"},{"containerName":null,"line":119,"kind":13,"localvar":"my","definition":"my","name":"$coords_q"},{"name":"$feat_r","kind":13,"line":127,"containerName":null},{"name":"$dbh","kind":13,"containerName":null,"line":127},{"kind":12,"line":127,"containerName":"main::","name":"selectcol_arrayref"},{"name":"$coords_q","kind":13,"containerName":null,"line":127},{"containerName":null,"line":131,"kind":13,"name":"%feat_r"},{"name":"$fid","definition":"my","line":132,"containerName":null,"kind":13,"localvar":"my"},{"line":132,"containerName":null,"kind":13,"name":"%feat_r"},{"kind":13,"containerName":null,"line":133,"name":"%DONE"},{"name":"%fid","containerName":null,"line":133,"kind":13},{"name":"$featr","definition":"my","containerName":null,"line":134,"kind":13,"localvar":"my"},{"name":"get_features_by_feature_id","line":134,"kind":12},{"name":"$dbh","kind":13,"containerName":null,"line":134},{"kind":13,"line":134,"containerName":null,"name":"$fid"},{"line":136,"containerName":null,"kind":13,"name":"%DONE"},{"kind":13,"containerName":null,"line":136,"name":"$fid"},{"kind":13,"containerName":null,"line":137,"name":"%featr"},{"kind":13,"containerName":null,"line":137,"name":"%fid"},{"name":"%seq_id","kind":13,"line":137,"containerName":null},{"name":"$min","containerName":null,"line":137,"kind":13},{"name":"%featr","containerName":null,"line":138,"kind":13},{"containerName":null,"line":138,"kind":13,"name":"%fid"},{"containerName":null,"line":138,"kind":13,"name":"%seq_id"},{"name":"$max","kind":13,"containerName":null,"line":138},{"containerName":null,"line":139,"kind":13,"name":"%featr"},{"kind":13,"line":139,"containerName":null,"name":"%fid"},{"line":139,"containerName":null,"kind":13,"name":"%seq_id"},{"name":"$strand","containerName":null,"line":139,"kind":13},{"name":"%featr","kind":13,"containerName":null,"line":140},{"name":"%fid","kind":13,"line":140,"containerName":null},{"kind":13,"line":140,"containerName":null,"name":"%seq_id"},{"name":"$min_partial","containerName":null,"line":140,"kind":13},{"line":141,"containerName":null,"kind":13,"name":"%featr"},{"line":141,"containerName":null,"kind":13,"name":"%fid"},{"kind":13,"line":141,"containerName":null,"name":"%seq_id"},{"kind":13,"containerName":null,"line":141,"name":"%max_partial"},{"containerName":null,"line":143,"kind":13,"name":"%featr"},{"name":"%fid","kind":13,"containerName":null,"line":143},{"name":"%seq_id","containerName":null,"line":143,"kind":13},{"kind":13,"line":144,"containerName":null,"name":"%featr"},{"name":"%fid","containerName":null,"line":144,"kind":13},{"kind":13,"line":144,"containerName":null,"name":"%seq_id"},{"name":"update_feature_mapping","line":146,"kind":12},{"containerName":null,"line":146,"kind":13,"name":"$dbh"},{"name":"$seq_id","containerName":null,"line":146,"kind":13},{"name":"%fid","kind":13,"containerName":null,"line":146},{"kind":13,"line":146,"containerName":null,"name":"$min"},{"name":"$max","kind":13,"line":147,"containerName":null},{"name":"$strand","kind":13,"line":148,"containerName":null},{"name":"$FT_coords","containerName":null,"line":149,"kind":13},{"name":"$min_partial","kind":13,"line":150,"containerName":null},{"containerName":null,"line":151,"kind":13,"name":"$max_partial"},{"containerName":null,"line":153,"kind":13,"name":"$feature"},{"name":"annotation","containerName":"main::","line":153,"kind":12},{"line":153,"containerName":"main::","kind":12,"name":"get_all_annotation_keys"},{"line":154,"kind":12,"name":"delete_feature_annotations"},{"name":"$dbh","line":154,"containerName":null,"kind":13},{"name":"%fid","kind":13,"containerName":null,"line":154},{"name":"$source","kind":13,"containerName":null,"line":154},{"line":155,"kind":12,"name":"load_feature_annotations"},{"line":155,"containerName":null,"kind":13,"name":"$dbh"},{"containerName":null,"line":155,"kind":13,"name":"$fid"},{"line":155,"containerName":null,"kind":13,"name":"$feature"},{"name":"annotation","containerName":"main::","line":155,"kind":12},{"name":"%source","kind":13,"line":155,"containerName":null},{"line":162,"containerName":null,"localvar":"my","kind":13,"name":"$feat_id","definition":"my"},{"line":162,"kind":12,"name":"load_SeqFeature"},{"kind":13,"line":162,"containerName":null,"name":"$dbh"},{"containerName":null,"line":162,"kind":13,"name":"$seq_id"},{"name":"$feature","line":162,"containerName":null,"kind":13},{"containerName":null,"line":162,"kind":13,"name":"$dbseqobj"},{"name":"$source","kind":13,"line":162,"containerName":null},{"name":"$locus_tag","kind":13,"line":163,"containerName":null},{"name":"$feature","containerName":null,"line":163,"kind":13},{"containerName":"main::","line":163,"kind":12,"name":"start"},{"line":163,"containerName":null,"kind":13,"name":"$feature"},{"kind":12,"containerName":"main::","line":163,"name":"end"},{"line":164,"containerName":null,"kind":13,"name":"$feature"},{"name":"annotation","kind":12,"line":164,"containerName":"main::"},{"kind":12,"containerName":"main::","line":164,"name":"get_all_annotation_keys"},{"line":165,"kind":12,"name":"load_feature_annotations"},{"containerName":null,"line":165,"kind":13,"name":"$dbh"},{"name":"$feat_id","containerName":null,"line":165,"kind":13},{"line":165,"containerName":null,"kind":13,"name":"$feature"},{"name":"annotation","containerName":"main::","line":165,"kind":12},{"name":"$source","containerName":null,"line":165,"kind":13},{"name":"%DONE","kind":13,"containerName":null,"line":167},{"kind":13,"containerName":null,"line":167,"name":"$feat_id"}]}