{"version":5,"vars":[{"name":"strict","kind":2,"line":2,"containerName":""},{"name":"GFF","line":3,"containerName":"Bio::Tools","kind":2},{"line":4,"kind":2,"containerName":"Bio::SeqFeature","name":"Generic"},{"line":5,"kind":2,"containerName":"Bio::DB","name":"Fasta"},{"name":"Std","containerName":"Getopt","line":6,"kind":2},{"name":"%ENV","kind":13,"line":7,"containerName":null},{"name":"ENVSCRIPTS","line":7,"kind":12},{"line":7,"kind":2,"containerName":"","name":"lib"},{"kind":2,"line":9,"containerName":"","name":"ENV"},{"kind":13,"definition":"my","line":11,"containerName":null,"name":"%opts","localvar":"my"},{"kind":12,"line":12,"name":"getopts"},{"name":"$opts","containerName":null,"line":12,"kind":13},{"name":"$dbh","line":16,"containerName":null,"localvar":"my","kind":13,"definition":"my"},{"name":"$annotation_file","containerName":null,"line":16,"kind":13},{"line":16,"containerName":null,"kind":13,"name":"$source"},{"kind":12,"line":16,"name":"handle_options"},{"name":"$opts","containerName":null,"line":16,"kind":13},{"definition":"my","kind":13,"localvar":"my","name":"$in","line":21,"containerName":null},{"line":21,"containerName":null,"kind":13,"name":"$annotation_file"},{"localvar":"my","line":22,"containerName":null,"name":"$line","definition":"my","kind":13},{"line":22,"containerName":null,"kind":13,"name":"%in"},{"name":"$line","line":23,"kind":13,"containerName":null},{"name":"$line","line":24,"containerName":null,"kind":13},{"kind":13,"definition":"my","name":"$query","line":25,"containerName":null,"localvar":"my"},{"name":"$seed_ortholog","kind":13,"line":26,"containerName":null},{"name":"$evalue","kind":13,"line":27,"containerName":null},{"containerName":null,"line":28,"kind":13,"name":"$score"},{"name":"$eggNOG_OGs","line":29,"kind":13,"containerName":null},{"containerName":null,"line":30,"kind":13,"name":"$max_annot_lvl"},{"name":"$COG_category","containerName":null,"line":31,"kind":13},{"kind":13,"line":32,"containerName":null,"name":"$description"},{"name":"$preferred_name","line":33,"containerName":null,"kind":13},{"name":"$GOs","containerName":null,"line":34,"kind":13},{"name":"$EC","line":35,"kind":13,"containerName":null},{"name":"$KEGG_ko","line":36,"kind":13,"containerName":null},{"kind":13,"line":37,"containerName":null,"name":"$KEGG_pathway"},{"containerName":null,"line":38,"kind":13,"name":"$KEGG_module"},{"name":"$KEGG_reaction","line":39,"kind":13,"containerName":null},{"name":"$KEGG_rclass","line":40,"containerName":null,"kind":13},{"kind":13,"line":41,"containerName":null,"name":"$BRITE"},{"name":"$KEGG_TC","containerName":null,"line":42,"kind":13},{"name":"$CAZy","kind":13,"line":43,"containerName":null},{"name":"$BiGG_reaction","containerName":null,"line":44,"kind":13},{"kind":13,"line":45,"containerName":null,"name":"$PFAMs"},{"line":45,"kind":13,"containerName":null,"name":"$line"},{"name":"@feat_r","line":49,"containerName":null,"localvar":"my","kind":13,"definition":"my"},{"kind":13,"line":50,"containerName":null,"name":"%query"},{"kind":13,"definition":"my","line":51,"containerName":null,"name":"$acc_q","localvar":"my"},{"kind":13,"line":52,"containerName":null,"name":"$feat_r"},{"name":"$dbh","line":52,"containerName":null,"kind":13},{"line":52,"kind":12,"containerName":"main::","name":"selectcol_arrayref"},{"name":"$acc_q","line":52,"containerName":null,"kind":13},{"kind":13,"line":56,"containerName":null,"name":"%feat_r"},{"definition":"my","kind":13,"localvar":"my","name":"$feat_id","line":57,"containerName":null},{"kind":13,"line":57,"containerName":null,"name":"@feat_r"},{"containerName":null,"line":58,"kind":13,"name":"%feat_r"},{"line":59,"kind":13,"containerName":null,"name":"$feat_r"},{"name":"delete_feature_annotations","line":62,"kind":12},{"line":62,"kind":13,"containerName":null,"name":"$dbh"},{"kind":13,"line":62,"containerName":null,"name":"%feat_id"},{"kind":12,"line":63,"name":"delete_feature_evidence"},{"name":"$dbh","containerName":null,"line":63,"kind":13},{"name":"%feat_id","line":63,"containerName":null,"kind":13},{"localvar":"my","name":"$feature_evidence_i","containerName":null,"line":64,"definition":"my","kind":13},{"kind":13,"definition":"my","name":"$feature_annotation_i","containerName":null,"line":68,"localvar":"my"},{"localvar":"my","name":"@NOGs","containerName":null,"line":75,"definition":"my","kind":13},{"line":75,"kind":13,"containerName":null,"name":"$eggNOG_OGs"},{"kind":13,"definition":"my","name":"%SEEN","containerName":null,"line":76,"localvar":"my"},{"definition":"my","kind":13,"localvar":"my","name":"$nog","line":77,"containerName":null},{"name":"@NOGs","containerName":null,"line":77,"kind":13},{"kind":13,"definition":"my","line":78,"containerName":null,"name":"@p","localvar":"my"},{"name":"$nog","line":78,"containerName":null,"kind":13},{"line":79,"kind":13,"containerName":null,"name":"@p"},{"name":"%SEEN","line":80,"kind":13,"containerName":null},{"name":"@p","kind":13,"line":80,"containerName":null},{"definition":"my","kind":13,"localvar":"my","name":"$prod","line":81,"containerName":null},{"containerName":null,"line":81,"kind":13,"name":"$dbh"},{"kind":12,"line":81,"containerName":"main::","name":"selectcol_arrayref"},{"name":"$description","line":82,"containerName":null,"kind":13},{"containerName":null,"line":82,"kind":13,"name":"@prod"},{"name":"$dbh","containerName":null,"line":83,"kind":13},{"containerName":"main::","line":83,"kind":12,"name":"do"},{"containerName":null,"line":83,"kind":13,"name":"%feature_evidence_i"},{"line":84,"containerName":null,"kind":13,"name":"@p"},{"line":84,"containerName":null,"kind":13,"name":"$DBI"},{"line":84,"containerName":"","kind":12,"name":"errstr"},{"name":"%SEEN","line":85,"kind":13,"containerName":null},{"name":"@p","kind":13,"line":85,"containerName":null},{"name":"%description","containerName":null,"line":91,"kind":13},{"name":"$dbh","kind":13,"line":92,"containerName":null},{"kind":12,"line":92,"containerName":"main::","name":"do"},{"line":92,"kind":13,"containerName":null,"name":"%feature_annotation_i"},{"line":93,"containerName":null,"kind":13,"name":"$description"},{"kind":13,"line":93,"containerName":null,"name":"$DBI"},{"name":"errstr","line":93,"kind":12,"containerName":""},{"name":"%preferred_name","line":96,"kind":13,"containerName":null},{"containerName":null,"line":97,"kind":13,"name":"$dbh"},{"name":"do","containerName":"main::","line":97,"kind":12},{"name":"%feature_annotation_i","line":97,"kind":13,"containerName":null},{"name":"$preferred_name","line":98,"kind":13,"containerName":null},{"name":"$DBI","line":98,"kind":13,"containerName":null},{"name":"errstr","line":98,"kind":12,"containerName":""},{"line":101,"containerName":null,"kind":13,"name":"%GOs"},{"kind":13,"definition":"my","name":"@gos","line":102,"containerName":null,"localvar":"my"},{"name":"$GOs","containerName":null,"line":102,"kind":13},{"name":"$go","line":103,"containerName":null,"localvar":"my","kind":13,"definition":"my"},{"containerName":null,"line":103,"kind":13,"name":"@gos"},{"kind":13,"line":104,"containerName":null,"name":"$dbh"},{"name":"do","line":104,"containerName":"main::","kind":12},{"name":"%feature_annotation_i","containerName":null,"line":104,"kind":13},{"name":"$go","containerName":null,"line":105,"kind":13},{"containerName":null,"line":105,"kind":13,"name":"$DBI"},{"name":"errstr","containerName":"","line":105,"kind":12},{"name":"%EC","line":109,"containerName":null,"kind":13},{"kind":13,"definition":"my","line":110,"containerName":null,"name":"@ecs","localvar":"my"},{"line":110,"kind":13,"containerName":null,"name":"$EC"},{"name":"$ec","line":111,"containerName":null,"localvar":"my","kind":13,"definition":"my"},{"line":111,"kind":13,"containerName":null,"name":"@ecs"},{"kind":13,"line":112,"containerName":null,"name":"$dbh"},{"line":112,"kind":12,"containerName":"main::","name":"do"},{"name":"%feature_annotation_i","containerName":null,"line":112,"kind":13},{"name":"$ec","line":113,"containerName":null,"kind":13},{"name":"$DBI","containerName":null,"line":113,"kind":13},{"containerName":"","line":113,"kind":12,"name":"errstr"},{"kind":13,"line":117,"containerName":null,"name":"%KEGG_ko"},{"kind":13,"definition":"my","containerName":null,"line":118,"name":"@kos","localvar":"my"},{"name":"$KEGG_ko","kind":13,"line":118,"containerName":null},{"definition":"my","kind":13,"localvar":"my","line":119,"containerName":null,"name":"$ko"},{"kind":13,"line":119,"containerName":null,"name":"@kos"},{"name":"$dbh","line":120,"containerName":null,"kind":13},{"name":"do","kind":12,"line":120,"containerName":"main::"},{"kind":13,"line":120,"containerName":null,"name":"%feature_evidence_i"},{"kind":13,"line":121,"containerName":null,"name":"$ko"},{"containerName":null,"line":121,"kind":13,"name":"$DBI"},{"containerName":"","line":121,"kind":12,"name":"errstr"},{"name":"%KEGG_TC","line":125,"containerName":null,"kind":13},{"localvar":"my","line":126,"containerName":null,"name":"@tcs","definition":"my","kind":13},{"name":"$KEGG_TC","line":126,"kind":13,"containerName":null},{"kind":13,"definition":"my","line":127,"containerName":null,"name":"$tc","localvar":"my"},{"line":127,"containerName":null,"kind":13,"name":"@tcs"},{"name":"$dbh","line":128,"containerName":null,"kind":13},{"name":"do","kind":12,"line":128,"containerName":"main::"},{"name":"%feature_annotation_i","line":128,"kind":13,"containerName":null},{"line":129,"kind":13,"containerName":null,"name":"$tc"},{"line":129,"kind":13,"containerName":null,"name":"$DBI"},{"containerName":"","line":129,"kind":12,"name":"errstr"},{"kind":13,"line":133,"containerName":null,"name":"%CAZy"},{"localvar":"my","name":"@cazys","line":134,"containerName":null,"definition":"my","kind":13},{"containerName":null,"line":134,"kind":13,"name":"$CAZy"},{"name":"$c","line":135,"containerName":null,"localvar":"my","kind":13,"definition":"my"},{"containerName":null,"line":135,"kind":13,"name":"@cazys"},{"line":136,"kind":13,"containerName":null,"name":"$dbh"},{"name":"do","containerName":"main::","line":136,"kind":12},{"containerName":null,"line":136,"kind":13,"name":"%feature_annotation_i"},{"name":"$c","line":137,"kind":13,"containerName":null},{"line":137,"kind":13,"containerName":null,"name":"$DBI"},{"name":"errstr","line":137,"containerName":"","kind":12},{"definition":"sub","kind":12,"children":[{"containerName":"handle_options","line":144,"name":"$opts","localvar":"my","kind":13,"definition":"my"},{"definition":"my","kind":13,"localvar":"my","name":"$err_msg","line":145,"containerName":"handle_options"},{"definition":"my","kind":13,"localvar":"my","name":"$dbh","line":147,"containerName":"handle_options"},{"name":"$opts","kind":13,"line":147,"containerName":"handle_options"},{"localvar":"my","line":148,"containerName":"handle_options","name":"$input_file","definition":"my","kind":13},{"name":"$opts","line":148,"containerName":"handle_options","kind":13},{"line":148,"kind":13,"containerName":"handle_options","name":"$err_msg"},{"localvar":"my","line":149,"containerName":"handle_options","name":"$source","definition":"my","kind":13},{"line":149,"containerName":"handle_options","kind":13,"name":"$opts"},{"containerName":"handle_options","line":151,"kind":13,"name":"$err_msg"},{"name":"$err_msg","line":151,"containerName":"handle_options","kind":13},{"name":"$dbh","containerName":"handle_options","line":153,"kind":13},{"name":"$input_file","line":153,"containerName":"handle_options","kind":13},{"kind":13,"line":153,"containerName":"handle_options","name":"$source"}],"range":{"end":{"line":154,"character":9999},"start":{"character":0,"line":143}},"name":"handle_options","line":143,"containerName":"main::"},{"name":"get_feat_types","containerName":"main::","line":156,"range":{"end":{"line":161,"character":9999},"start":{"line":156,"character":0}},"kind":12,"children":[{"localvar":"my","line":157,"containerName":"get_feat_types","name":"$dbh","definition":"my","kind":13},{"kind":13,"definition":"my","name":"$q","line":158,"containerName":"get_feat_types","localvar":"my"},{"kind":13,"definition":"my","line":159,"containerName":"get_feat_types","name":"$r","localvar":"my"},{"line":159,"containerName":"get_feat_types","kind":13,"name":"$dbh"},{"line":159,"kind":12,"containerName":"get_feat_types","name":"selectall_hashref"},{"kind":13,"line":159,"containerName":"get_feat_types","name":"$q"},{"containerName":"get_feat_types","line":160,"kind":13,"name":"$r"}],"definition":"sub"}]}