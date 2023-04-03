# Some convenient definitions of selectors that we can use in our configuration.
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos

selections = {
    't->Wb': parsel.MultiSelection( # t->Wb from ttbar
        [
            parsel.FirstSelector(22,6), # top quark
            parsel.FirstSelector(23,5), # bottom quark
            parsel.FirstSelector(22,24) # W boson
        ]
    ),
    'Wb': parsel.MultiSelection( # Wb from ttbar
        [
            parsel.FirstSelector(23,5), # bottom quark
            parsel.FirstSelector(22,24) # W boson
        ]
    ),
    't->Wb w/ qq and W daughters': parsel.MultiSelection(
        [
            parsel.FirstSelector(22, 6), # top quark
            parsel.FirstSelector(23, 5), # bottom quark
            parsel.FirstSelector(22,24), # W boson
            parsel.AlgoSelection(algos.SelectSimplestQuarks(     parsel.FirstSelector(22,24)),n=2, fixed_length=True), # q's from W->qq'
            parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(22,24)),n=120) # up to 120 stable daughters of W
        ]
    ),
    't daughters': parsel.MultiSelection(
        [
            parsel.AlgoSelection(algos.SelectFinalStateDaughters(parsel.FirstSelector(22,6)),n=200) # up to 200 stable daughters of top quark
        ]
    ),

}