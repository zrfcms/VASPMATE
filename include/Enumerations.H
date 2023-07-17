
/** \name Font Numbers
    The following constants define the standard FLTK fonts:
 */
/*@{*/
/** A font number is an index into the internal font table. */
typedef int Fl_Font;

const Fl_Font FL_HELVETICA              = 0;	///< Helvetica (or Arial) normal (0)
const Fl_Font FL_HELVETICA_BOLD         = 1;	///< Helvetica (or Arial) bold
const Fl_Font FL_HELVETICA_ITALIC       = 2;	///< Helvetica (or Arial) oblique
const Fl_Font FL_HELVETICA_BOLD_ITALIC  = 3;	///< Helvetica (or Arial) bold-oblique
const Fl_Font FL_COURIER                = 4;	///< Courier normal
const Fl_Font FL_COURIER_BOLD           = 5;	///< Courier bold 
const Fl_Font FL_COURIER_ITALIC         = 6;	///< Courier italic
const Fl_Font FL_COURIER_BOLD_ITALIC    = 7;	///< Courier bold-italic
const Fl_Font FL_TIMES                  = 8;	///< Times roman
const Fl_Font FL_TIMES_BOLD             = 9;	///< Times roman bold
const Fl_Font FL_TIMES_ITALIC           = 10;	///< Times roman italic
const Fl_Font FL_TIMES_BOLD_ITALIC      = 11;	///< Times roman bold-italic
const Fl_Font FL_SYMBOL                 = 12;	///< Standard symbol font
const Fl_Font FL_SCREEN                 = 13;	///< Default monospaced screen font
const Fl_Font FL_SCREEN_BOLD            = 14;	///< Default monospaced bold screen font
const Fl_Font FL_ZAPF_DINGBATS          = 15;	///< Zapf-dingbats font

const Fl_Font FL_FREE_FONT              = 16;	///< first one to allocate
const Fl_Font FL_BOLD                   = 1;	///< add this to helvetica, courier, or times
const Fl_Font FL_ITALIC                 = 2;	///< add this to helvetica, courier, or times
const Fl_Font FL_BOLD_ITALIC            = 3;	///< add this to helvetica, courier, or times

/*@}*/

/** Size of a font in pixels.
    This is the approximate height of a font in pixels.
 */
typedef int Fl_Fontsize;

/** \name Colors 
    The Fl_Color type holds an FLTK color value.

    Colors are either 8-bit indexes into a <a href="fltk-colormap.png">virtual colormap</a>
    or 24-bit RGB color values. (See \ref drawing_colors for the default FLTK colormap)

    Color indices occupy the lower 8 bits of the value, while
    RGB colors occupy the upper 24 bits, for a byte organization of RGBI.

<pre>
 Fl_Color => 0xrrggbbii
                | | | |
                | | | +--- \ref drawing_colors "index" between 0 and 255
                | | +----- blue color component (8 bit)
                | +------- green component (8 bit)
                +--------- red component (8 bit)
 </pre>

    A color can have either an index or an rgb value. Colors with rgb set 
    and an index >0 are reserved for special use.

 */
/*@{*/
/** An FLTK color value; see also \ref drawing_colors  */
typedef unsigned int Fl_Color;

// Standard colors. These are used as default colors in widgets and altered as necessary
const Fl_Color FL_FOREGROUND_COLOR  = 0;	///< the default foreground color (0) used for labels and text
const Fl_Color FL_BACKGROUND2_COLOR = 7;	///< the default background color for text, list, and valuator widgets
const Fl_Color FL_INACTIVE_COLOR    = 8;	///< the inactive foreground color
const Fl_Color FL_SELECTION_COLOR   = 15;	///< the default selection/highlight color

  // boxtypes generally limit themselves to these colors so
  // the whole ramp is not allocated:

const Fl_Color FL_GRAY0   = 32;			// 'A'
const Fl_Color FL_DARK3   = 39;			// 'H'
const Fl_Color FL_DARK2   = 45;			// 'N'
const Fl_Color FL_DARK1   = 47;			// 'P'
const Fl_Color FL_BACKGROUND_COLOR  = 49;	// 'R' default background color
const Fl_Color FL_LIGHT1  = 50;			// 'S'
const Fl_Color FL_LIGHT2  = 52;			// 'U'
const Fl_Color FL_LIGHT3  = 54;			// 'W'

  // FLTK provides a 5x8x5 color cube that is used with colormap visuals

const Fl_Color FL_BLACK   = 56;
const Fl_Color FL_RED     = 88;
const Fl_Color FL_GREEN   = 63;
const Fl_Color FL_YELLOW  = 95;
const Fl_Color FL_BLUE    = 216;
const Fl_Color FL_MAGENTA = 248;
const Fl_Color FL_CYAN    = 223;
const Fl_Color FL_DARK_RED = 72;

const Fl_Color FL_DARK_GREEN    = 60;
const Fl_Color FL_DARK_YELLOW   = 76;
const Fl_Color FL_DARK_BLUE     = 136;
const Fl_Color FL_DARK_MAGENTA  = 152;
const Fl_Color FL_DARK_CYAN     = 140;

const Fl_Color FL_WHITE         = 255;