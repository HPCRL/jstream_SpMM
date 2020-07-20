#ifdef PAPI
	int event_set = PAPI_NULL;
	long_long values[4];
	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	values [0] = values[1] = values [2] = values[3] = 0;
	
	if(retval != PAPI_VER_CURRENT && retval > 0){
		cout<<retval << " "<< PAPI_VER_CURRENT<<endl;
	}

	string event_names[4];

	event_names[0] = event_names[1] = event_names[2] = event_names[3] = "event";


	if ((retval = PAPI_create_eventset(&event_set)) != PAPI_OK) {
		fprintf(stderr, "PAPI error %d: %s\n",retval,           PAPI_strerror(retval));
		exit(1);
	}

	//PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_L1_DCM));
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_L2_DCM));
	event_names[0] += "PAPI_L2_DCM";
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_L3_TCM));
	event_names[1] += "PAPI_L3_TCM";
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_PRF_DM));
	event_names[2] += "PAPI_PRF_DM";

	PAPI_ERROR_CHECK(PAPI_start(event_set));
#endif