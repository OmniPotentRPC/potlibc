/* SPDX-FileCopyrightText: 2026 Rohit Goswami */
/* SPDX-License-Identifier: BSD-3-Clause */

/* Bit goldens of the CBE gradient jet, in Bohr and Hartree.
 *
 * Generated from the first-order-to-third-order forward-mode kernel of the
 * spectral-instanton ac071977 reference build (libcbe.so cbe_eval, cbe_hvp and
 * cbe_jet), reproduced bit for bit by potlibc's ch52008 surface at the saddle point and at a displaced geometry. Every value is
 * the exact IEEE-754 double the kernel returns. gradient is dE/dx (the force
 * is its negative); hvp, cubic and quartic are the first three derivatives
 * with respect to t of the gradient along the stated direction. */
/* saddle */
static const double saddle_x[18] = {
    0x1.4102d30bb0c6ep+1,    0x1.85a1650848229p+1,
    0x1.4d4e7ba318765p+1,    0x1.1c2f62bef5f7dp+0,
    0x1.5c81ba817f3dap+0,    0x1.260b8652234aap+0,
    0x1.c6d9ace433becp-7,    0x1.ab480071c99ecp-3,
    0x1.3a8fe3cce0444p+1,    0x1.3b07a8b639e16p+1,
    0x1.e30cc3be763cep-3,    0x1.3546d4ad7b93ap-4,
    -0x1.8be26c89dfe6dp-4,    0x1.3f31a9f2a98dp+1,
    -0x1.4091239c1dafbp-4,    0x1.b1534bc5f0994p+1,
    0x1.066a688e7452bp+2,    0x1.c23e9d64ad619p+1,
};
static const double saddle_v[18] = {
    0x1.e921dd42f09bap-1,    -0x1.de67ac55f1633p-7,
    -0x1.f1216dba340c9p-2,    -0x1.f606eec8ac71dp-2,
    0x1.6ad6c3c07d448p-2,    0x1.bd21aaf88dcfbp-2,
    -0x1.f2c746fb6af5cp-3,    -0x1.ffd7c43328f13p-2,
    -0x1.2a47bccaac8bep-3,    0x1.b00da046b65ebp-1,
    0x1.7c499e9bc960cp-2,    -0x1.c93355db185ecp-3,
    -0x1.f696931df044ap-1,    -0x1.424411148142bp-5,
    0x1.e1096828da6fdp-2,    0x1.29a3033463914p-1,
    -0x1.41cd40f14f49fp-2,    -0x1.d5ccca4d1fa0cp-2,
};
static const double saddle_energy = 0x1.87b76b3d8e2cp-6;
static const double saddle_gradient[18] = {
    -0x1.77652a4cec941p-39,    -0x1.0bb24b61fa787p-37,
    -0x1.cbc5b7e6f92cbp-41,    0x1.efb7940a99af7p-38,
    0x1.78edd730160adp-36,    -0x1.a7a47c46f4d1dp-38,
    0x1.960683b6ba1dap-40,    -0x1.095492a2a0557p-40,
    0x1.dfd3d2721ee82p-38,    0x1.e30e5a13ef7f8p-40,
    -0x1.0e66634a72a39p-38,    0x1.33556f59858f2p-39,
    -0x1.9c2a22c12e573p-38,    -0x1.43fcd37b27bc1p-38,
    -0x1.92daa867f0a87p-39,    -0x1.d87edfb523097p-40,
    -0x1.379ab1b2babadp-38,    0x1.8860effc15d9fp-41,
};
static const double saddle_hvp[18] = {
    0x1.7f97887638bfdp-5,    -0x1.bdede93cb7feap-8,
    -0x1.aae607052bbe5p-6,    -0x1.12984fa8b11dcp-2,
    0x1.019c9337881e5p-2,    0x1.136f59c6f3901p-2,
    -0x1.0289171bbb056p-6,    -0x1.022bd8d4cd823p-7,
    -0x1.c8f3e59ad7f61p-6,    0x1.2b28418e2a58ep-2,
    -0x1.aeadc764d56efp-3,    -0x1.88bc2f1322e1p-3,
    -0x1.c70377871b2fap-6,    -0x1.1d21dc85b2a77p-9,
    -0x1.920df9269992cp-8,    -0x1.bea1a0a12f005p-6,
    -0x1.902556078b57fp-6,    -0x1.18b6b8ec77dd8p-6,
};
static const double saddle_cubic[18] = {
    -0x1.98d6c840026dbp-3,    -0x1.7baceebdb040cp-4,
    -0x1.cfd63931e87acp-5,    0x1.09ce1a3efbc17p+0,
    -0x1.ceff23df50d03p-1,    -0x1.86bd0cac197eep-1,
    -0x1.80ffaf6deafd3p-4,    -0x1.e37a6a93293a4p-5,
    0x1.c1bb26d7c62c6p-4,    -0x1.7c45ff6ba20dcp-1,
    0x1.f1d810d5b6af7p-1,    0x1.5620a2ac1995bp-1,
    -0x1.08ac4428e249ep-4,    0x1.5a4b3c4082db2p-7,
    -0x1.6cee41283a723p-7,    0x1.00a7db8426743p-4,
    0x1.2b5954cc058b6p-4,    0x1.b1621bcc694c6p-5,
};
static const double saddle_quartic[18] = {
    0x1.f9669edc40b87p-7,    0x1.3c236e810c8fdp-2,
    0x1.c61489d1744a5p-2,    -0x1.4de88462aa3a1p-2,
    0x1.9055730374c5dp+1,    0x1.825ebad9fa7d8p-1,
    0x1.c119835acd74dp-4,    -0x1.6e3851a138f0bp-3,
    0x1.6e5e7e6485e39p-5,    0x1.92fa143919443p-3,
    -0x1.89455bf5681f2p+1,    -0x1.26632da4312f4p+0,
    0x1.fb4b13441f74fp-4,    -0x1.53fd410cdfb6bp-5,
    -0x1.b003779a9759fp-8,    -0x1.e9e38161fea27p-4,
    -0x1.2610abfe72aeap-3,    -0x1.61442b8b2b966p-4,
};
static const long saddle_hvp_operations = 17172;
static const double saddle_hvp_scale = 0x1.2b28418e2a58ep-2;
static const long saddle_jet_operations = 17172;
static const double saddle_jet_scale = 0x1.9055730374c5dp+1;
/* displaced */
static const double displaced_x[18] = {
    0x1.46657d8eb2964p+1,    0x1.8bfa236684dabp+1,
    0x1.51a1296dd939dp+1,    0x1.1cb7a32239c5p+0,
    0x1.54accae7ed57p+0,    0x1.19885a9fd3974p+0,
    -0x1.f04d083a1d036p-6,    0x1.84ff0c9e83a35p-3,
    0x1.3c8e51b8112bap+1,    0x1.407903af0c924p+1,
    0x1.242e1945e2d0dp-2,    0x1.bd0d830ea606ep-4,
    -0x1.86cf7a6838a61p-4,    0x1.3b318c50cb921p+1,
    -0x1.04bb396918c99p-3,    0x1.abb8f48473ce5p+1,
    0x1.0544f26e4ef52p+2,    0x1.c4572620d62ap+1,
};
static const double displaced_v[18] = {
    0x1.e921dd42f09bap-1,    -0x1.de67ac55f1633p-7,
    -0x1.f1216dba340c9p-2,    -0x1.f606eec8ac71dp-2,
    0x1.6ad6c3c07d448p-2,    0x1.bd21aaf88dcfbp-2,
    -0x1.f2c746fb6af5cp-3,    -0x1.ffd7c43328f13p-2,
    -0x1.2a47bccaac8bep-3,    0x1.b00da046b65ebp-1,
    0x1.7c499e9bc960cp-2,    -0x1.c93355db185ecp-3,
    -0x1.f696931df044ap-1,    -0x1.424411148142bp-5,
    0x1.e1096828da6fdp-2,    0x1.29a3033463914p-1,
    -0x1.41cd40f14f49fp-2,    -0x1.d5ccca4d1fa0cp-2,
};
static const double displaced_energy = 0x1.967c4c9a341cp-6;
static const double displaced_gradient[18] = {
    -0x1.fe5c76ed7c98fp-10,    -0x1.0626c11433c17p-8,
    -0x1.bc8967547ba92p-9,    0x1.e40d57179355cp-6,
    0x1.c4efa2d794601p-10,    -0x1.31e63961d89e8p-6,
    -0x1.7cf2825ac34a5p-7,    -0x1.1029618720112p-7,
    0x1.a48d7f68e308dp-7,    -0x1.a8e8d2a45ee28p-7,
    0x1.af3d5f7b5a4d1p-7,    0x1.5ba9b6345fc78p-7,
    -0x1.e3751d79cf374p-11,    -0x1.86c7035996ecep-14,
    0x1.4e1339d47d35bp-11,    -0x1.10f1e1eadfe03p-9,
    -0x1.46440ef97ee4bp-9,    -0x1.08a67286e9e9cp-9,
};
static const double displaced_hvp[18] = {
    0x1.9426a5228a1efp-5,    0x1.052ebe8dc8783p-6,
    0x1.2d1d3864cd444p-11,    -0x1.588eb2fc3c748p-2,
    0x1.2eec0618a190dp-2,    0x1.33f91aef3ae9bp-2,
    -0x1.065f5b6b0283ep-7,    -0x1.5b6ebe2a1954bp-7,
    -0x1.13f93a9c694b7p-5,    0x1.70b3c31a6ef4fp-2,
    -0x1.07925e9bd4d94p-2,    -0x1.d310319dc3fc9p-3,
    -0x1.a051a2209a943p-6,    -0x1.327bfa4ae0169p-9,
    -0x1.5d070555c1572p-8,    -0x1.438e7e2910356p-5,
    -0x1.53612bfe15a36p-5,    -0x1.18a26a9d3919ap-5,
};
static const double displaced_cubic[18] = {
    -0x1.3957a3abace6dp-3,    -0x1.cdd4c62a5d851p-5,
    -0x1.14a1f2d280ce1p-6,    0x1.35af043e96642p+0,
    -0x1.20128e8c5fe36p+0,    -0x1.b2cc83c64b8e7p-1,
    -0x1.350542feb0d28p-4,    -0x1.49fb4b193f2e1p-5,
    0x1.7361433aacbe1p-4,    -0x1.ea3a1d099d956p-1,
    0x1.28b813e2d135ap+0,    0x1.8b547b39bb36ep-1,
    -0x1.14db18db781d4p-4,    0x1.26bf1b2e0745p-6,
    -0x1.e4762f4a0d0dfp-7,    0x1.66e08f2a1247p-5,
    0x1.6fbfd8de6ec54p-5,    0x1.28590f1edf52fp-6,
};
static const double displaced_quartic[18] = {
    0x1.cfbfd1cace8c2p-5,    0x1.3ed13fe51f2a7p-2,
    0x1.baac4c75d109bp-2,    -0x1.8bb7d085dc47ap-1,
    0x1.ef2c05d909bbap+1,    0x1.06a492ada3c94p+0,
    0x1.181300d0674ecp-4,    -0x1.23019c6c087bcp-3,
    0x1.0a2f4bac4d6adp-4,    0x1.4507ff956abc4p-1,
    -0x1.ea369bbacd5d7p+1,    -0x1.66ba7ad6e3b68p+0,
    0x1.df13681415f31p-4,    -0x1.f669c85d6bb5ep-5,
    0x1.2212f492e45e8p-8,    -0x1.a987ca46582e6p-4,
    -0x1.2c5d132aa0d13p-3,    -0x1.02d1951c60838p-3,
};
static const long displaced_hvp_operations = 17172;
static const double displaced_hvp_scale = 0x1.70b3c31a6ef4fp-2;
static const long displaced_jet_operations = 17172;
static const double displaced_jet_scale = 0x1.ef2c05d909bbap+1;
