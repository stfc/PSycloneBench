import "regent"

terralib.includepath = ".;"..terralib.includepath
local stdlib = terralib.includec("stdlib.h")
local stdio = terralib.includec("stdio.h")
local gocean2d_io_mod = terralib.includec("gocean2d_io_mod.h")
terralib.linklibrary("libgfortran.so")
terralib.linklibrary("./libgocean2d_io_mod.so")

fspace setup_type{
  jpiglo: int, 
  jpjglo: int, 
  dx: double,
  dy: double,
  nit000: int,
  nitend: int,
  record: int,
  jphgr_msh: int,
  dep_const: double,
  rdt: double,
  cbfr: double,
  visc: double
}

--Uses the libgocean2d_io_mod to read the namelist into the setup_type
task read_namelist(in1 : region(ispace(int1d), setup_type)) where writes(in1) do
  var t_jpi = 0
  var t_jpj = 0
  var t_dx : double = 0.0 
  var t_dy : double = 0.0
  var t_nit000 = 0
  var t_nitend = 0
  var t_record = 0
  var t_jphgr = 0
  var t_dep : double = 0.0
  var t_rdt : double = 0.0
  var t_cbfr : double = 0.0
  var t_visc : double = 0.0
  gocean2d_io_mod.read_namelist(&t_jpi, &t_jpj, &t_dx, &t_dy, &t_nit000, &t_nitend,
                                &t_record, &t_jphgr, &t_dep, &t_rdt, &t_cbfr, &t_visc)
  in1[0].jpiglo = t_jpi
  in1[0].jpjglo = t_jpj
  in1[0].dx = t_dx
  in1[0].dy = t_dy
  in1[0].nit000 = t_nit000
  in1[0].nitend = t_nitend
  in1[0].record = t_record
  in1[0].jphgr_msh = t_jphgr
  in1[0].dep_const = t_dep
  in1[0].rdt = t_rdt
  in1[0].cbfr = t_cbfr
  in1[0].visc = t_visc
end
 
task print_setup_data(output: region(ispace(int1d), setup_type))
where reads(output) do
 stdio.printf("%i %i %f %f %i %i %i %i %f %f %f %f\n", output[0].jpiglo, output[0].jpjglo, output[0].dx,
              output[0].dy, output[0].nit000, output[0].nitend, output[0].record, output[0].jphgr_msh,
              output[0].dep_const, output[0].rdt, output[0].cbfr, output[0].visc); 
end


task main()

  var setup_data_is = ispace(int1d, 1)
  var setup_data = region(setup_data_is, setup_type)
  fill(setup_data.{jpiglo, jpjglo, nit000, nitend, record, jphgr_msh}, 0)
  fill(setup_data.{dx, dy, dep_const, rdt, cbfr, visc}, 0)
--  init(setup_data)

--  print(setup_data)

end

--regentlib.start(main)