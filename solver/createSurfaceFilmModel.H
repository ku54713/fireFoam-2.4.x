Info<< "\nConstructing surface film model" << nl;

typedef regionModels::surfaceFilmModels::surfaceFilmModel filmModelType;

autoPtr<filmModelType> tsurfaceFilm(filmModelType::New(mesh, g));
filmModelType& surfaceFilm = tsurfaceFilm();
