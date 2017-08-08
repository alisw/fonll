if (!defined &_FPU_CONTROL_H) {
    eval 'sub _FPU_CONTROL_H () {1;}' unless defined(&_FPU_CONTROL_H);
    require 'features.ph';
    eval 'sub _FPU_MASK_IM () {0x01;}' unless defined(&_FPU_MASK_IM);
    eval 'sub _FPU_MASK_DM () {0x02;}' unless defined(&_FPU_MASK_DM);
    eval 'sub _FPU_MASK_ZM () {0x04;}' unless defined(&_FPU_MASK_ZM);
    eval 'sub _FPU_MASK_OM () {0x08;}' unless defined(&_FPU_MASK_OM);
    eval 'sub _FPU_MASK_UM () {0x10;}' unless defined(&_FPU_MASK_UM);
    eval 'sub _FPU_MASK_PM () {0x20;}' unless defined(&_FPU_MASK_PM);
    eval 'sub _FPU_EXTENDED () {0x300;}' unless defined(&_FPU_EXTENDED);
    eval 'sub _FPU_DOUBLE () {0x200;}' unless defined(&_FPU_DOUBLE);
    eval 'sub _FPU_SINGLE () {0x0;}' unless defined(&_FPU_SINGLE);
    eval 'sub _FPU_RC_NEAREST () {0x0;}' unless defined(&_FPU_RC_NEAREST);
    eval 'sub _FPU_RC_DOWN () {0x400;}' unless defined(&_FPU_RC_DOWN);
    eval 'sub _FPU_RC_UP () {0x800;}' unless defined(&_FPU_RC_UP);
    eval 'sub _FPU_RC_ZERO () {0xC00;}' unless defined(&_FPU_RC_ZERO);
    eval 'sub _FPU_RESERVED () {0xF0C0;}' unless defined(&_FPU_RESERVED);
    eval 'sub _FPU_DEFAULT () {0x137f;}' unless defined(&_FPU_DEFAULT);
    eval 'sub _FPU_IEEE () {0x137f;}' unless defined(&_FPU_IEEE);
    eval 'sub _FPU_GETCW {
        local($cw) = @_;
        eval " &__asm__ (\\"fnstcw %0\\" : \\"=m\\" (*$cw))";
    }' unless defined(&_FPU_GETCW);
    eval 'sub _FPU_SETCW {
        local($cw) = @_;
        eval " &__asm__ (\\"fldcw %0\\" : : \\"m\\" (*$cw))";
    }' unless defined(&_FPU_SETCW);
}
1;
