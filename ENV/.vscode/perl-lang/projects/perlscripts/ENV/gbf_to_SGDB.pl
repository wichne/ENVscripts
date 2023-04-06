{"vars":[{"kind":2,"line":2,"containerName":"","name":"DBI"},{"kind":2,"containerName":"Bio","line":3,"name":"SeqIO"},{"line":4,"containerName":"Bio","name":"Seq","kind":2},{"kind":2,"containerName":"Getopt","line":5,"name":"Std"},{"name":"Dumper","containerName":"Data","line":6,"kind":2},{"name":"strict","containerName":"","line":7,"kind":2},{"line":8,"containerName":null,"name":"%ENV","kind":13},{"line":8,"name":"SCRIPTS","kind":12},{"kind":2,"containerName":"","line":8,"name":"lib"},{"kind":2,"name":"SGDB","line":10,"containerName":""},{"kind":2,"name":"vars","line":11,"containerName":""},{"name":"%args","containerName":null,"kind":13,"definition":"my","line":13,"localvar":"my"},{"name":"getopts","line":14,"kind":12},{"kind":13,"containerName":null,"line":14,"name":"$args"},{"name":"$filename","containerName":null,"kind":13,"definition":"my","localvar":"my","line":16},{"kind":13,"name":"%args","containerName":null,"line":16},{"kind":13,"line":17,"localvar":"my","definition":"my","containerName":null,"name":"$format"},{"name":"%args","containerName":null,"line":17,"kind":13},{"kind":13,"line":17,"containerName":null,"name":"%args"},{"containerName":null,"name":"$db","line":18,"localvar":"my","definition":"my","kind":13},{"name":"%args","line":18,"containerName":null,"kind":13},{"kind":13,"definition":"my","localvar":"my","line":19,"name":"$pswd","containerName":null},{"name":"%args","line":19,"containerName":null,"kind":13},{"kind":13,"definition":"my","localvar":"my","line":20,"name":"$user","containerName":null},{"containerName":null,"line":20,"name":"%args","kind":13},{"name":"%args","line":20,"containerName":null,"kind":13},{"containerName":null,"line":20,"name":"%ENV","kind":13},{"name":"$asmbl_id","containerName":null,"kind":13,"definition":"my","line":21,"localvar":"my"},{"kind":13,"name":"%args","line":21,"containerName":null},{"containerName":null,"name":"$in","line":23,"localvar":"my","definition":"my","kind":13},{"name":"Bio","containerName":"SeqIO","line":23,"kind":12},{"name":"new","line":23,"containerName":"main::","kind":12},{"kind":13,"line":23,"containerName":null,"name":"$filename"},{"kind":13,"name":"$format","line":24,"containerName":null},{"definition":"my","line":26,"localvar":"my","kind":13,"name":"$host","containerName":null},{"kind":13,"containerName":null,"line":26,"name":"%ENV"},{"name":"DBSERVER","line":26,"kind":12},{"name":"%ENV","line":26,"containerName":null,"kind":13},{"kind":12,"line":26,"name":"DBSERVER"},{"definition":"my","localvar":"my","line":27,"kind":13,"name":"$user","containerName":null},{"line":27,"containerName":null,"name":"%ENV","kind":13},{"kind":12,"name":"USER","line":27},{"kind":13,"definition":"my","line":28,"localvar":"my","name":"$dbh","containerName":null},{"kind":12,"line":28,"name":"DBI"},{"containerName":"main::","line":28,"name":"connect","kind":12},{"containerName":null,"line":28,"name":"$user","kind":13},{"name":"%pswd","containerName":null,"line":28,"kind":13},{"kind":13,"name":"%args","line":30,"containerName":null},{"name":"@tables","containerName":null,"kind":13,"definition":"my","localvar":"my","line":31},{"name":"$t","containerName":null,"definition":"my","line":38,"localvar":"my","kind":13},{"containerName":null,"line":38,"name":"@tables","kind":13},{"name":"$q","containerName":null,"kind":13,"definition":"my","localvar":"my","line":39},{"kind":13,"line":40,"containerName":null,"name":"$dbh"},{"kind":12,"name":"do","containerName":"main::","line":40},{"containerName":null,"line":40,"name":"$q","kind":13},{"definition":"my","localvar":"my","line":44,"kind":13,"name":"$q","containerName":null},{"line":45,"localvar":"my","definition":"my","kind":13,"containerName":null,"name":"$existing"},{"name":"$dbh","line":45,"containerName":null,"kind":13},{"line":45,"containerName":"main::","name":"selectall_hashref","kind":12},{"name":"$q","containerName":null,"line":45,"kind":13},{"kind":13,"line":47,"localvar":"my","definition":"my","containerName":null,"name":"$seqo"},{"kind":13,"name":"$in","containerName":null,"line":47},{"kind":12,"name":"next_seq","line":47,"containerName":"main::"},{"definition":"my","localvar":"my","line":48,"kind":13,"name":"%FEAT","containerName":null},{"definition":"my","line":49,"localvar":"my","kind":13,"name":"@features","containerName":null},{"line":49,"containerName":null,"name":"$seqo","kind":13},{"kind":12,"name":"get_SeqFeatures","containerName":"main::","line":49},{"kind":13,"definition":"my","line":52,"localvar":"my","name":"$seqf","containerName":null},{"kind":13,"line":52,"containerName":null,"name":"@features"},{"kind":13,"containerName":null,"line":53,"name":"$asmbl_id"},{"name":"load_assembly","line":53,"kind":12},{"kind":13,"line":53,"containerName":null,"name":"$seqf"},{"kind":13,"name":"%args","containerName":null,"line":53},{"containerName":null,"line":54,"name":"$seqo","kind":13},{"kind":12,"name":"display_id","line":54,"containerName":"main::"},{"kind":13,"line":55,"localvar":"my","definition":"my","containerName":null,"name":"$featObj"},{"kind":13,"containerName":null,"line":55,"name":"@features"},{"kind":12,"line":56,"name":"process_feature"},{"containerName":null,"line":56,"name":"$asmbl_id","kind":13},{"name":"$featObj","line":56,"containerName":null,"kind":13},{"line":56,"containerName":null,"name":"%FEAT","kind":13},{"name":"$locus","containerName":null,"kind":13,"definition":"my","line":58,"localvar":"my"},{"kind":13,"line":58,"containerName":null,"name":"%FEAT"},{"kind":13,"line":59,"containerName":null,"name":"%existing"},{"kind":13,"line":59,"containerName":null,"name":"$locus"},{"name":"$feat_type","containerName":null,"kind":13,"definition":"my","line":60,"localvar":"my"},{"name":"%FEAT","containerName":null,"line":60,"kind":13},{"line":60,"containerName":null,"name":"%locus","kind":13},{"kind":13,"definition":"my","line":61,"localvar":"my","name":"$feat_name","containerName":null},{"kind":12,"name":"get_next_feat_name","line":61},{"kind":13,"name":"%dbh","containerName":null,"line":61},{"containerName":null,"line":61,"name":"$feat_type","kind":13},{"kind":13,"line":62,"containerName":null,"name":"%FEAT"},{"containerName":null,"line":62,"name":"%locus","kind":13},{"name":"$feat_name","containerName":null,"line":62,"kind":13},{"kind":12,"line":63,"name":"insert_feature"},{"line":63,"containerName":null,"name":"$dbh","kind":13},{"containerName":null,"line":63,"name":"%FEAT","kind":13},{"name":"%locus","containerName":null,"line":63,"kind":13},{"name":"%FEAT","containerName":null,"line":64,"kind":13},{"kind":13,"line":64,"containerName":null,"name":"%locus"},{"kind":13,"name":"$feat_name","containerName":null,"line":64},{"kind":12,"name":"insert_ident","line":65},{"kind":13,"name":"$dbh","containerName":null,"line":65},{"line":65,"containerName":null,"name":"%FEAT","kind":13},{"kind":13,"name":"%locus","line":65,"containerName":null},{"kind":13,"containerName":null,"line":66,"name":"%FEAT"},{"containerName":null,"line":66,"name":"%locus","kind":13},{"kind":13,"name":"$feat_name","line":66,"containerName":null},{"kind":12,"line":67,"name":"insert_translation"},{"kind":13,"line":67,"containerName":null,"name":"$dbh"},{"containerName":null,"line":67,"name":"%FEAT","kind":13},{"kind":13,"name":"%locus","containerName":null,"line":67},{"name":"load_assembly","containerName":"main::","range":{"end":{"line":96,"character":9999},"start":{"character":0,"line":76}},"definition":"sub","line":76,"kind":12,"children":[{"definition":"my","line":77,"localvar":"my","kind":13,"name":"$assemblyObj","containerName":"load_assembly"},{"containerName":"load_assembly","name":"%assembly","localvar":"my","line":79,"definition":"my","kind":13},{"kind":13,"containerName":"load_assembly","line":79,"name":"$assemblyObj"},{"kind":12,"name":"seq","line":79,"containerName":"load_assembly"},{"kind":12,"line":79,"containerName":"load_assembly","name":"seq"},{"name":"$asmbl_id","containerName":"load_assembly","kind":13,"definition":"my","line":80,"localvar":"my"},{"kind":13,"containerName":"load_assembly","line":80,"name":"$dbh"},{"name":"%assembly","line":80,"containerName":"load_assembly","kind":13},{"name":"%asmbl_data","containerName":"load_assembly","definition":"my","line":82,"localvar":"my","kind":13},{"kind":13,"line":82,"containerName":"load_assembly","name":"$assemblyObj"},{"containerName":"load_assembly","line":82,"name":"seq","kind":12},{"line":82,"containerName":"load_assembly","name":"description","kind":12},{"kind":13,"line":83,"containerName":"load_assembly","name":"$assemblyObj"},{"kind":12,"line":83,"containerName":"load_assembly","name":"get_tag_values"},{"line":87,"localvar":"my","definition":"my","kind":13,"containerName":"load_assembly","name":"$asmbl_data_id"},{"name":"$dbh","containerName":"load_assembly","line":87,"kind":13},{"name":"%asmbl_data","containerName":"load_assembly","line":87,"kind":13},{"containerName":"load_assembly","name":"%stan","kind":13,"localvar":"my","line":89,"definition":"my"},{"kind":13,"containerName":"load_assembly","line":89,"name":"$asmbl_id"},{"name":"$asmbl_data_id","containerName":"load_assembly","line":90,"kind":13},{"containerName":"load_assembly","line":91,"name":"$asmbl_id","kind":13},{"name":"$dbh","line":94,"containerName":"load_assembly","kind":13},{"kind":13,"name":"%stan","line":94,"containerName":"load_assembly"},{"line":95,"containerName":"load_assembly","name":"$asmbl_id","kind":13}]},{"kind":12,"line":80,"name":"insert_assembly"},{"kind":12,"line":87,"name":"insert_asmbl_data"},{"name":"insert_stan","line":94,"kind":12},{"range":{"start":{"line":98,"character":0},"end":{"character":9999,"line":243}},"containerName":"main::","name":"process_feature","kind":12,"children":[{"line":102,"localvar":"my","definition":"my","kind":13,"containerName":"process_feature","name":"$asmbl_id"},{"kind":13,"definition":"my","line":103,"localvar":"my","name":"$featObj","containerName":"process_feature"},{"containerName":"process_feature","name":"$FEAT","kind":13,"localvar":"my","line":104,"definition":"my"},{"definition":"my","line":106,"localvar":"my","kind":13,"name":"$feat_type","containerName":"process_feature"},{"kind":13,"line":106,"containerName":"process_feature","name":"$featObj"},{"kind":12,"name":"primary_tag","line":106,"containerName":"process_feature"},{"name":"@tags","containerName":"process_feature","kind":13,"definition":"my","line":107,"localvar":"my"},{"kind":13,"containerName":"process_feature","line":107,"name":"$featObj"},{"kind":12,"name":"get_all_tags","containerName":"process_feature","line":107},{"kind":13,"name":"@tags","line":164,"containerName":"process_feature"},{"kind":13,"name":"$feat_type","containerName":"process_feature","line":168},{"name":"$locus","containerName":"process_feature","definition":"my","localvar":"my","line":169,"kind":13},{"containerName":"process_feature","line":169,"name":"$featObj","kind":13},{"kind":12,"name":"get_tag_values","containerName":"process_feature","line":169},{"containerName":"process_feature","name":"$end5","localvar":"my","line":170,"definition":"my","kind":13},{"kind":13,"line":170,"containerName":"process_feature","name":"$end3"},{"containerName":"process_feature","line":170,"name":"$featObj","kind":13},{"containerName":"process_feature","line":170,"name":"strand","kind":12},{"kind":13,"name":"$featObj","containerName":"process_feature","line":171},{"containerName":"process_feature","line":171,"name":"start","kind":12},{"line":171,"containerName":"process_feature","name":"$featObj","kind":13},{"name":"end","containerName":"process_feature","line":171,"kind":12},{"line":172,"containerName":"process_feature","name":"$featObj","kind":13},{"containerName":"process_feature","line":172,"name":"end","kind":12},{"name":"$featObj","containerName":"process_feature","line":172,"kind":13},{"kind":12,"line":172,"containerName":"process_feature","name":"start"},{"line":174,"containerName":"process_feature","name":"$feat_type","kind":13},{"name":"$feat_type","containerName":"process_feature","line":174,"kind":13},{"name":"$FEAT","containerName":"process_feature","line":175,"kind":13},{"kind":13,"name":"$locus","containerName":"process_feature","line":175},{"containerName":"process_feature","line":175,"name":"$asmbl_id","kind":13},{"line":176,"containerName":"process_feature","name":"$FEAT","kind":13},{"kind":13,"line":176,"containerName":"process_feature","name":"$locus"},{"kind":13,"line":176,"containerName":"process_feature","name":"$feat_type"},{"containerName":"process_feature","line":177,"name":"$FEAT","kind":13},{"kind":13,"containerName":"process_feature","line":177,"name":"$locus"},{"kind":13,"containerName":"process_feature","line":177,"name":"$end5"},{"name":"$FEAT","containerName":"process_feature","line":178,"kind":13},{"kind":13,"containerName":"process_feature","line":178,"name":"$locus"},{"name":"$end3","containerName":"process_feature","line":178,"kind":13},{"kind":13,"name":"$FEAT","containerName":"process_feature","line":179},{"line":179,"containerName":"process_feature","name":"$locus","kind":13},{"name":"$featObj","line":179,"containerName":"process_feature","kind":13},{"containerName":"process_feature","line":179,"name":"seq","kind":12},{"kind":12,"name":"seq","containerName":"process_feature","line":179},{"kind":13,"containerName":"process_feature","line":180,"name":"$FEAT"},{"kind":13,"name":"$locus","containerName":"process_feature","line":180},{"kind":13,"name":"$featObj","containerName":"process_feature","line":180},{"name":"source_tag","line":180,"containerName":"process_feature","kind":12},{"kind":13,"line":182,"localvar":"my","definition":"my","containerName":"process_feature","name":"$tag"},{"name":"@tags","line":182,"containerName":"process_feature","kind":13},{"containerName":"process_feature","name":"@values","localvar":"my","line":183,"definition":"my","kind":13},{"name":"$featObj","containerName":"process_feature","line":183,"kind":13},{"kind":12,"name":"get_tag_values","containerName":"process_feature","line":183},{"line":183,"containerName":"process_feature","name":"$tag","kind":13},{"kind":13,"containerName":"process_feature","line":184,"name":"$tag"},{"kind":13,"containerName":"process_feature","line":184,"name":"$FEAT"},{"kind":13,"name":"$locus","containerName":"process_feature","line":184},{"kind":13,"line":184,"containerName":"process_feature","name":"@values"},{"name":"$tag","containerName":"process_feature","line":185,"kind":13},{"kind":13,"name":"$FEAT","containerName":"process_feature","line":185},{"kind":13,"name":"$locus","containerName":"process_feature","line":185},{"line":185,"containerName":"process_feature","name":"@values","kind":13},{"name":"$tag","line":186,"containerName":"process_feature","kind":13},{"containerName":"process_feature","line":186,"name":"$FEAT","kind":13},{"kind":13,"containerName":"process_feature","line":186,"name":"$locus"},{"name":"$values","containerName":"process_feature","line":186,"kind":13},{"kind":13,"name":"$tag","line":187,"containerName":"process_feature"},{"line":187,"containerName":"process_feature","name":"$FEAT","kind":13},{"name":"$locus","containerName":"process_feature","line":187,"kind":13},{"kind":13,"name":"@values","line":187,"containerName":"process_feature"},{"kind":13,"name":"$tag","containerName":"process_feature","line":188},{"name":"$FEAT","containerName":"process_feature","line":195,"kind":13},{"line":195,"containerName":"process_feature","name":"$locus","kind":13},{"kind":13,"containerName":"process_feature","line":195,"name":"$values"},{"kind":13,"containerName":"process_feature","line":197,"name":"$tag"},{"kind":13,"name":"$FEAT","containerName":"process_feature","line":198},{"kind":13,"name":"$locus","line":198,"containerName":"process_feature"},{"kind":13,"containerName":"process_feature","line":198,"name":"$values"},{"containerName":"process_feature","line":199,"name":"$FEAT","kind":13},{"containerName":"process_feature","line":199,"name":"$locus","kind":13},{"line":199,"containerName":"process_feature","name":"$values","kind":13},{"name":"$tag","containerName":"process_feature","line":200,"kind":13},{"kind":13,"name":"$FEAT","containerName":"process_feature","line":200},{"line":200,"containerName":"process_feature","name":"$locus","kind":13},{"name":"$values","containerName":"process_feature","line":200,"kind":13},{"name":"$tag","line":201,"containerName":"process_feature","kind":13},{"localvar":"my","line":202,"definition":"my","kind":13,"containerName":"process_feature","name":"$v"},{"containerName":"process_feature","line":202,"name":"@values","kind":13},{"kind":13,"containerName":"process_feature","line":203,"name":"$FEAT"},{"containerName":"process_feature","line":203,"name":"$locus","kind":13},{"kind":13,"name":"$tag","line":206,"containerName":"process_feature"},{"kind":13,"line":206,"containerName":"process_feature","name":"$FEAT"},{"name":"$locus","line":206,"containerName":"process_feature","kind":13},{"name":"$tag","containerName":"process_feature","line":207,"kind":13},{"line":207,"containerName":"process_feature","name":"$FEAT","kind":13},{"kind":13,"name":"$locus","line":207,"containerName":"process_feature"},{"kind":13,"line":208,"containerName":"process_feature","name":"$tag"},{"name":"$FEAT","containerName":"process_feature","line":208,"kind":13},{"name":"$locus","containerName":"process_feature","line":208,"kind":13},{"containerName":"process_feature","line":208,"name":"@values","kind":13},{"kind":13,"name":"$tag","containerName":"process_feature","line":209},{"kind":13,"name":"$FEAT","containerName":"process_feature","line":209},{"line":209,"containerName":"process_feature","name":"$locus","kind":13},{"kind":13,"containerName":"process_feature","line":212,"name":"$FEAT"},{"kind":13,"containerName":"process_feature","line":212,"name":"$locus"},{"kind":13,"name":"$end5","containerName":"process_feature","line":212},{"kind":13,"containerName":"process_feature","line":212,"name":"$end3"},{"containerName":"process_feature","line":212,"name":"$FEAT","kind":13},{"kind":13,"name":"$locus","line":212,"containerName":"process_feature"},{"kind":13,"containerName":"process_feature","line":212,"name":"$FEAT"},{"name":"$locus","line":212,"containerName":"process_feature","kind":13},{"kind":13,"name":"$FEAT","line":212,"containerName":"process_feature"},{"containerName":"process_feature","line":212,"name":"$locus","kind":13},{"line":214,"containerName":"process_feature","name":"$featObj","kind":13},{"kind":12,"name":"primary_tag","line":214,"containerName":"process_feature"},{"kind":13,"definition":"my","localvar":"my","line":215,"name":"@coords","containerName":"process_feature"},{"kind":13,"containerName":"process_feature","line":216,"name":"$featObj"},{"name":"location","containerName":"process_feature","line":216,"kind":12},{"containerName":"process_feature","line":216,"name":"isa","kind":12},{"kind":13,"definition":"my","line":217,"localvar":"my","name":"$location","containerName":"process_feature"},{"line":217,"containerName":"process_feature","name":"$featObj","kind":13},{"containerName":"process_feature","line":217,"name":"location","kind":12},{"kind":12,"containerName":"process_feature","line":217,"name":"sub_Location"},{"line":218,"containerName":"process_feature","name":"$location","kind":13},{"kind":12,"name":"strand","line":218,"containerName":"process_feature"},{"name":"@coords","containerName":"process_feature","line":219,"kind":13},{"containerName":"process_feature","line":219,"name":"$location","kind":13},{"name":"start","line":219,"containerName":"process_feature","kind":12},{"containerName":"process_feature","line":219,"name":"$location","kind":13},{"kind":12,"name":"end","containerName":"process_feature","line":219},{"containerName":"process_feature","line":221,"name":"@coords","kind":13},{"kind":13,"line":221,"containerName":"process_feature","name":"$location"},{"name":"end","containerName":"process_feature","line":221,"kind":12},{"name":"$location","line":221,"containerName":"process_feature","kind":13},{"kind":12,"name":"start","containerName":"process_feature","line":221},{"name":"@coords","containerName":"process_feature","line":225,"kind":13},{"name":"$featObj","line":225,"containerName":"process_feature","kind":13},{"line":225,"containerName":"process_feature","name":"strand","kind":12},{"containerName":"process_feature","line":226,"name":"$featObj","kind":13},{"kind":12,"line":226,"containerName":"process_feature","name":"start"},{"name":"$featObj","line":226,"containerName":"process_feature","kind":13},{"kind":12,"line":226,"containerName":"process_feature","name":"end"},{"name":"$featObj","containerName":"process_feature","line":227,"kind":13},{"containerName":"process_feature","line":227,"name":"end","kind":12},{"kind":13,"name":"$featObj","containerName":"process_feature","line":227},{"kind":12,"name":"start","line":227,"containerName":"process_feature"},{"line":229,"containerName":"process_feature","name":"$FEAT","kind":13},{"kind":13,"line":229,"containerName":"process_feature","name":"$locus"},{"kind":13,"name":"@coords","containerName":"process_feature","line":229},{"name":"$FEAT","line":232,"containerName":"process_feature","kind":13},{"line":232,"containerName":"process_feature","name":"$locus","kind":13},{"containerName":"process_feature","line":234,"name":"$FEAT","kind":13},{"kind":13,"name":"$locus","containerName":"process_feature","line":234},{"name":"$FEAT","line":236,"containerName":"process_feature","kind":13},{"kind":13,"line":236,"containerName":"process_feature","name":"$locus"},{"name":"$FEAT","line":236,"containerName":"process_feature","kind":13},{"kind":13,"name":"$locus","containerName":"process_feature","line":236},{"line":239,"containerName":"process_feature","name":"$FEAT","kind":13},{"name":"$locus","line":239,"containerName":"process_feature","kind":13}],"line":98,"definition":"sub"},{"kind":12,"name":"next","line":168},{"range":{"end":{"line":312,"character":9999},"start":{"character":0,"line":245}},"name":"load_feature","containerName":"main::","children":[{"name":"$asmbl_id","containerName":"load_feature","kind":13,"definition":"my","localvar":"my","line":246},{"kind":13,"line":247,"localvar":"my","definition":"my","containerName":"load_feature","name":"$featObj"},{"name":"$ident","containerName":"load_feature","kind":13,"definition":"my","line":249,"localvar":"my"},{"line":249,"containerName":"load_feature","name":"$transl","kind":13},{"name":"$featObj","containerName":"load_feature","line":249,"kind":13},{"containerName":"load_feature","name":"$feat_type","kind":13,"line":251,"localvar":"my","definition":"my"},{"kind":13,"name":"$featObj","line":251,"containerName":"load_feature"},{"kind":12,"line":251,"containerName":"load_feature","name":"primary_tag"},{"containerName":"load_feature","name":"$feat_name","line":253,"localvar":"my","definition":"my","kind":13},{"kind":13,"line":253,"containerName":"load_feature","name":"$feat_type"},{"name":"$end5","containerName":"load_feature","kind":13,"definition":"my","line":256,"localvar":"my"},{"kind":13,"containerName":"load_feature","line":256,"name":"$end3"},{"line":256,"containerName":"load_feature","name":"$featObj","kind":13},{"kind":12,"line":256,"containerName":"load_feature","name":"strand"},{"kind":13,"name":"$featObj","containerName":"load_feature","line":257},{"name":"start","containerName":"load_feature","line":257,"kind":12},{"name":"$featObj","line":257,"containerName":"load_feature","kind":13},{"name":"end","line":257,"containerName":"load_feature","kind":12},{"name":"$featObj","containerName":"load_feature","line":258,"kind":13},{"line":258,"containerName":"load_feature","name":"end","kind":12},{"name":"$featObj","containerName":"load_feature","line":258,"kind":13},{"kind":12,"containerName":"load_feature","line":258,"name":"start"},{"name":"$ident","line":260,"containerName":"load_feature","kind":13},{"name":"%feature","containerName":"load_feature","kind":13,"definition":"my","localvar":"my","line":264},{"containerName":"load_feature","line":264,"name":"$asmbl_id","kind":13},{"kind":13,"containerName":"load_feature","line":265,"name":"$feat_type"},{"containerName":"load_feature","line":266,"name":"$feat_name","kind":13},{"containerName":"load_feature","line":267,"name":"$end5","kind":13},{"kind":13,"name":"$end3","line":268,"containerName":"load_feature"},{"name":"$featObj","line":269,"containerName":"load_feature","kind":13},{"kind":12,"containerName":"load_feature","line":269,"name":"seq"},{"kind":12,"name":"seq","containerName":"load_feature","line":269},{"kind":13,"name":"$featObj","containerName":"load_feature","line":270},{"line":270,"containerName":"load_feature","name":"source_tag","kind":12},{"name":"$feature","containerName":"load_feature","line":272,"kind":13},{"kind":13,"line":272,"containerName":"load_feature","name":"$transl"},{"kind":13,"name":"$transl","line":272,"containerName":"load_feature"},{"name":"$dbh","containerName":"load_feature","line":274,"kind":13},{"kind":13,"name":"%feature","line":274,"containerName":"load_feature"},{"containerName":"load_feature","line":276,"name":"$feat_type","kind":13},{"kind":13,"name":"$feat_type","line":278,"containerName":"load_feature"},{"kind":13,"containerName":"load_feature","line":279,"name":"$transl"},{"kind":13,"line":279,"containerName":"load_feature","name":"$feat_name"},{"name":"@coords","containerName":"load_feature","definition":"my","line":280,"localvar":"my","kind":13},{"containerName":"load_feature","line":281,"name":"$featObj","kind":13},{"line":281,"containerName":"load_feature","name":"location","kind":12},{"kind":12,"line":281,"containerName":"load_feature","name":"isa"},{"containerName":"load_feature","name":"$location","kind":13,"line":282,"localvar":"my","definition":"my"},{"name":"$featObj","containerName":"load_feature","line":282,"kind":13},{"containerName":"load_feature","line":282,"name":"location","kind":12},{"kind":12,"line":282,"containerName":"load_feature","name":"sub_Location"},{"name":"$location","containerName":"load_feature","line":283,"kind":13},{"name":"strand","containerName":"load_feature","line":283,"kind":12},{"kind":13,"containerName":"load_feature","line":284,"name":"@coords"},{"line":284,"containerName":"load_feature","name":"$location","kind":13},{"name":"start","line":284,"containerName":"load_feature","kind":12},{"kind":13,"containerName":"load_feature","line":284,"name":"$location"},{"kind":12,"containerName":"load_feature","line":284,"name":"end"},{"kind":13,"line":286,"containerName":"load_feature","name":"@coords"},{"containerName":"load_feature","line":286,"name":"$location","kind":13},{"containerName":"load_feature","line":286,"name":"end","kind":12},{"kind":13,"name":"$location","line":286,"containerName":"load_feature"},{"name":"start","containerName":"load_feature","line":286,"kind":12},{"name":"@coords","line":290,"containerName":"load_feature","kind":13},{"line":290,"containerName":"load_feature","name":"$featObj","kind":13},{"containerName":"load_feature","line":290,"name":"strand","kind":12},{"kind":13,"line":291,"containerName":"load_feature","name":"$featObj"},{"kind":12,"name":"start","containerName":"load_feature","line":291},{"name":"$featObj","containerName":"load_feature","line":291,"kind":13},{"kind":12,"line":291,"containerName":"load_feature","name":"end"},{"line":292,"containerName":"load_feature","name":"$featObj","kind":13},{"kind":12,"name":"end","containerName":"load_feature","line":292},{"kind":13,"name":"$featObj","line":292,"containerName":"load_feature"},{"kind":12,"containerName":"load_feature","line":292,"name":"start"},{"name":"$transl","line":294,"containerName":"load_feature","kind":13},{"kind":13,"name":"@coords","containerName":"load_feature","line":294},{"containerName":"load_feature","line":295,"name":"$transl","kind":13},{"name":"$transl","containerName":"load_feature","line":295,"kind":13},{"containerName":"load_feature","line":296,"name":"$dbh","kind":13},{"name":"$transl","containerName":"load_feature","line":296,"kind":13},{"kind":13,"name":"$ident","containerName":"load_feature","line":299},{"name":"$ident","containerName":"load_feature","line":301,"kind":13},{"kind":13,"containerName":"load_feature","line":303,"name":"$ident"},{"kind":13,"containerName":"load_feature","line":303,"name":"$ident"},{"name":"$ident","line":306,"containerName":"load_feature","kind":13},{"name":"$ident","line":309,"containerName":"load_feature","kind":13},{"name":"$feat_name","line":309,"containerName":"load_feature","kind":13},{"name":"$dbh","line":310,"containerName":"load_feature","kind":13},{"kind":13,"name":"$ident","line":310,"containerName":"load_feature"}],"kind":12,"definition":"sub","line":245},{"kind":12,"line":249,"name":"parse_tags"},{"kind":12,"name":"get_next_feat_name","line":253},{"name":"insert_feature","line":274,"kind":12},{"line":276,"name":"return","kind":12},{"kind":12,"name":"insert_translation","line":296},{"line":310,"name":"insert_ident","kind":12},{"name":"parse_tags","containerName":"main::","range":{"start":{"line":314,"character":0},"end":{"character":9999,"line":376}},"definition":"sub","line":314,"children":[{"containerName":"parse_tags","name":"$featObj","line":315,"localvar":"my","definition":"my","kind":13},{"name":"@tags","containerName":"parse_tags","definition":"my","localvar":"my","line":316,"kind":13},{"name":"$featObj","line":316,"containerName":"parse_tags","kind":13},{"containerName":"parse_tags","line":316,"name":"get_all_tags","kind":12},{"name":"%ident","containerName":"parse_tags","kind":13,"definition":"my","localvar":"my","line":318},{"name":"%transl","containerName":"parse_tags","definition":"my","line":319,"localvar":"my","kind":13},{"containerName":"parse_tags","name":"$locus","line":321,"localvar":"my","definition":"my","kind":13},{"name":"$protein","line":321,"containerName":"parse_tags","kind":13},{"name":"$product","containerName":"parse_tags","line":321,"kind":13},{"line":321,"containerName":"parse_tags","name":"$gene_sym","kind":13},{"containerName":"parse_tags","line":321,"name":"$ec","kind":13},{"kind":13,"name":"@GO","containerName":"parse_tags","line":321},{"kind":13,"name":"$comment","containerName":"parse_tags","line":321},{"localvar":"my","line":323,"definition":"my","kind":13,"containerName":"parse_tags","name":"$tag"},{"kind":13,"containerName":"parse_tags","line":323,"name":"@tags"},{"name":"@values","containerName":"parse_tags","kind":13,"definition":"my","localvar":"my","line":324},{"kind":13,"name":"$featObj","containerName":"parse_tags","line":324},{"kind":12,"line":324,"containerName":"parse_tags","name":"get_tag_values"},{"kind":13,"containerName":"parse_tags","line":324,"name":"$tag"},{"containerName":"parse_tags","line":326,"name":"$tag","kind":13},{"line":327,"containerName":"parse_tags","name":"$ident","kind":13},{"kind":13,"containerName":"parse_tags","line":327,"name":"@values"},{"line":330,"containerName":"parse_tags","name":"$tag","kind":13},{"kind":13,"containerName":"parse_tags","line":331,"name":"$ident"},{"kind":13,"line":331,"containerName":"parse_tags","name":"@values"},{"kind":13,"name":"$tag","line":334,"containerName":"parse_tags"},{"line":335,"containerName":"parse_tags","name":"$ident","kind":13},{"kind":13,"name":"$values","containerName":"parse_tags","line":335},{"name":"$tag","containerName":"parse_tags","line":338,"kind":13},{"name":"$ident","line":339,"containerName":"parse_tags","kind":13},{"containerName":"parse_tags","line":339,"name":"@values","kind":13},{"kind":13,"name":"$tag","containerName":"parse_tags","line":342},{"kind":13,"line":349,"containerName":"parse_tags","name":"$ident"},{"kind":13,"containerName":"parse_tags","line":349,"name":"$values"},{"line":352,"containerName":"parse_tags","name":"$tag","kind":13},{"kind":13,"name":"$transl","containerName":"parse_tags","line":353},{"kind":13,"line":353,"containerName":"parse_tags","name":"$values"},{"line":356,"containerName":"parse_tags","name":"$tag","kind":13},{"name":"$transl","line":356,"containerName":"parse_tags","kind":13},{"kind":13,"containerName":"parse_tags","line":356,"name":"$values"},{"line":358,"containerName":"parse_tags","name":"$tag","kind":13},{"definition":"my","localvar":"my","line":359,"kind":13,"name":"$v","containerName":"parse_tags"},{"kind":13,"name":"@values","containerName":"parse_tags","line":359},{"line":360,"containerName":"parse_tags","name":"$transl","kind":13},{"kind":13,"containerName":"parse_tags","line":364,"name":"$tag"},{"kind":13,"name":"$transl","line":364,"containerName":"parse_tags"},{"line":366,"containerName":"parse_tags","name":"$tag","kind":13},{"kind":13,"name":"$transl","line":366,"containerName":"parse_tags"},{"name":"$tag","containerName":"parse_tags","line":368,"kind":13},{"kind":13,"line":369,"containerName":"parse_tags","name":"$transl"},{"containerName":"parse_tags","line":369,"name":"@values","kind":13},{"name":"$tag","line":372,"containerName":"parse_tags","kind":13},{"kind":13,"name":"$transl","line":372,"containerName":"parse_tags"},{"kind":13,"name":"$ec","containerName":"parse_tags","line":374},{"kind":13,"name":"$ec","containerName":"parse_tags","line":374},{"kind":13,"name":"%ident","line":375,"containerName":"parse_tags"},{"kind":13,"line":375,"containerName":"parse_tags","name":"%transl"}],"kind":12}],"version":5}