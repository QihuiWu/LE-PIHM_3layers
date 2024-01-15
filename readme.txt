pihm.h
102   	int				bedrock;  /* bedrock type */

read_alloc.c
149   			fscanf(att_file, "%d %d %d %d", &(DS->Ele[i].soil), &(DS->Ele[i].geol), &(DS->Ele[i].bedrock), &(DS->Ele[i].LC));

Initialize.c
209			DS->Ele[i].RhoBed = DS->LE_bedrock[(DS->Ele[i].bedrock - 1)].RhoBed;
211  			DS->Ele[i].Uplift = DS->LE_bedrock[(DS->Ele[i].bedrock - 1)].Uplift; /*Bedrock uplift rate*/
212  			DS->Ele[i].CoefP0 = CS->Cal.CoefP0*DS->LE_bedrock[(DS->Ele[i].bedrock - 1)].CoefP0; /*Fitting constant in weather equation*/
214  			DS->Ele[i].P0 = CS->Cal.DReg*DS->LE_bedrock[(DS->Ele[i].bedrock - 1)].P0; /*Regolith production rate in the absence of soil above bedrock*/



