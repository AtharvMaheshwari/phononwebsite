export class LocalEventDB {
    constructor() {
        this.name = "Model Materials";
        this.year = 2026;
        this.author = "L. Wirtz";
        this.url = "#";
    }

    isAvailable() {
        return false;
    }

    get_materials(callback) {
        let reference = this.author+", "+"<a href='"+this.url+"'>"+this.name+"</a> ("+this.year+")";
        let name = this.name;

        function dothings(materials) {
            for (let i=0; i<materials.length; i++) {
                let m = materials[i];
                m.source = name;
                m.type = "json";
                m.reference = reference;
                m.url = "data/jsonfiles_event_LudgerWirtz/" + m.file;
            }
            callback(materials)
        }

        $.get('data/jsonfiles_event_LudgerWirtz/models.json', dothings);
    }
}
