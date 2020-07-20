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
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_RES_STL));
	event_names[0] += "PAPI_RES_STL";
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_MEM_WCY));
	event_names[1] += "PAPI_MEM_WCY";
	PAPI_ERROR_CHECK(PAPI_add_event(event_set, PAPI_FUL_ICY));
	event_names[2] += "PAPI_FUL_ICY";

	PAPI_ERROR_CHECK(PAPI_start(event_set));
#endif